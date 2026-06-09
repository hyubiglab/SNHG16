import pandas as pd
import pysam
import gffutils
import sys
from Bio.Seq import Seq
from collections import defaultdict
import pyfaidx
import re

P_SITE_OFFSETS = {}

# ORF codon settings (non canonical included)
START_CODONS = ['ATG', 'GTG', 'CTG', 'TTG']
STOP_CODONS = ['TGA', 'TAA', 'TAG']

# sORF settings
MIN_ORF_LENGTH_BP = 30
MAX_ORF_LENGTH_BP = 303

MIN_P_SITES_FOR_CODING = 3  # Minimum threashold for coding potential
MIN_FRAME_0_RATIO = 50.0

def get_gene_id_from_name(db, gene_name):
    gene_ids = []

    try:
        for gene in db.features_of_type('gene'):
            attrs = gene.attributes

            if 'gene_name' in attrs and attrs['gene_name'][0] == gene_name:
                gene_ids.append(gene.id)
            elif 'Name' in attrs and attrs['Name'][0] == gene_name:
                gene_ids.append(gene.id)
            elif gene.id == gene_name:
                gene_ids.append(gene.id)

    except Exception as e:
        print(f"[WARNING] gene_name search failed for {gene_name}: {e}")

    return gene_ids


def get_gene_name_from_gene_feature(gene, fallback):
    attrs = gene.attributes
    if 'gene_name' in attrs:
        return attrs['gene_name'][0]
    if 'Name' in attrs:
        return attrs['Name'][0]
    if 'gene_id' in attrs:
        return attrs['gene_id'][0]
    return fallback


def load_transcript_tpm_table(tpm_path):
    df = pd.read_csv(tpm_path, sep=r"\s+", engine="python")

    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce").fillna(0)
    df["percentage"] = pd.to_numeric(df["percentage"], errors="coerce").fillna(0)

    return df[["gene_name", "gene_id", "transcript_id", "TPM", "percentage"]]

def get_transcript_features(db, transcript_id):
    try:
        transcript = db[transcript_id]
        strand = transcript.strand
        seqid = transcript.seqid
    except gffutils.NoFeatError:
        return None, None, None

    exon_regions = list(db.children(transcript, featuretype='exon'))
    exon_regions.sort(key=lambda x: x.start)

    return exon_regions, strand, seqid


def load_transcript_sequence_from_genome(transcript_id, db, genome_fasta_path):
    exon_regions, strand, seqid = get_transcript_features(db, transcript_id)
    if not exon_regions:
        return None

    try:
        genome = pyfaidx.Fasta(genome_fasta_path)
    except Exception as e:
        print(f"[ERROR] failed to load genome FASTA: {e}")
        return None

    transcript_seq = []

    for exon in exon_regions:
        try:
            exon_seq = genome[exon.seqid][exon.start - 1:exon.end].seq
            transcript_seq.append(exon_seq)
        except KeyError:
            print(f"[WARNING] chromosome not found in FASTA: {exon.seqid}")
            return None

    final_seq = ''.join(transcript_seq).upper()

    if strand == '-':
        final_seq = str(Seq(final_seq).reverse_complement()).upper()

    return final_seq


def transcript_to_genomic_coordinate(exon_regions, ts_coord, strand):
    current_ts_len = 0
    exons = exon_regions if strand == '+' else exon_regions[::-1]

    for exon in exons:
        exon_len = exon.end - exon.start + 1

        if current_ts_len + exon_len > ts_coord:
            offset_in_exon = ts_coord - current_ts_len

            if strand == '+':
                return exon.start + offset_in_exon
            else:
                return exon.end - offset_in_exon

        current_ts_len += exon_len

    return None


def find_all_orf_candidates(seq, start_codons, stop_codons, min_length_bp):
    seq = seq.upper()
    all_orf_candidates = []

    for frame in range(3):
        for orf_start in range(frame, len(seq) - 2, 3):
            if seq[orf_start:orf_start + 3] not in start_codons:
                continue

            # Closest stop codon within same frame
            for stop_pos in range(orf_start + 3, len(seq) - 2, 3):
                if seq[stop_pos:stop_pos + 3] in stop_codons:
                    orf_end = stop_pos + 3
                    orf_length = orf_end - orf_start

                    if min_length_bp <= orf_length < MAX_ORF_LENGTH_BP:
                        all_orf_candidates.append((orf_start, orf_end))

                    break

    return sorted(list(set(all_orf_candidates)))


def calculate_psite_coverage(bam_file, seqid, start, end, strand, exon_regions=None):
    
    coverage = defaultdict(int)

    exon_positions = set()
    if exon_regions:
        for exon in exon_regions:
            exon_positions.update(range(exon.start, exon.end + 1))

    try:
        read_iter = bam_file.fetch(seqid, start, end)
    except ValueError as e:
        print(f"[WARNING] BAM fetch failed for {seqid}:{start}-{end}: {e}")
        return {}

    for read in read_iter:
        if read.is_unmapped:
            continue

        read_len = read.query_length
        if read_len not in P_SITE_OFFSETS:
            continue

        offset = P_SITE_OFFSETS[read_len]

        if read.is_reverse:
            p_site = read.reference_end - 1 - offset
        else:
            p_site = read.reference_start + offset

        if read.is_reverse != (strand == '-'):
            continue

        if exon_regions and p_site not in exon_positions:
            continue

        coverage[p_site] += 1

    return dict(coverage)


def calculate_periodicity_p_sites(orf_ts_start, orf_ts_end, bam_file, transcript_id, db):
    exon_regions, strand, seqid = get_transcript_features(db, transcript_id)
    if not exon_regions:
        return 0, 0.0, 0, [0, 0, 0]

    orf_genomic_start = transcript_to_genomic_coordinate(exon_regions, orf_ts_start, strand)
    orf_genomic_end = transcript_to_genomic_coordinate(exon_regions, orf_ts_end - 1, strand)

    if orf_genomic_start is None or orf_genomic_end is None:
        return 0, 0.0, 0, [0, 0, 0]

    g_start = min(e.start for e in exon_regions)
    g_end = max(e.end for e in exon_regions)

    genomic_coverage = calculate_psite_coverage(
        bam_file,
        seqid,
        g_start,
        g_end,
        strand,
        exon_regions=exon_regions
    )

    frame_counts = [0, 0, 0]

    for g_loci, count in genomic_coverage.items():
        if g_loci < min(orf_genomic_start, orf_genomic_end):
            continue
        if g_loci > max(orf_genomic_start, orf_genomic_end):
            continue

        if strand == '+':
            distance = g_loci - orf_genomic_start
        else:
            distance = orf_genomic_start - g_loci

        if distance < 0:
            continue

        frame = distance % 3
        frame_counts[frame] += count

    total_orf_p_sites = sum(frame_counts)
    frame_0_count = frame_counts[0]

    if total_orf_p_sites == 0:
        return 0, 0.0, 0, frame_counts

    frame_0_ratio = (frame_0_count / total_orf_p_sites) * 100

    return total_orf_p_sites, frame_0_ratio, frame_0_count, frame_counts


def choose_main_orf(orf_results):
    if len(orf_results) == 0:
        return None

    return max(
        orf_results,
        key=lambda x: (
            x.get("predicted_orf_p_sites", 0),
            x.get("frame_0_ratio (%)", 0),
            x.get("orf_length_bp", 0),
            -x.get("orf_ts_start", 10**12)
        )
    )


def classify_relative_orf(orf_start, orf_end, main_start, main_end):

    if orf_start == main_start and orf_end == main_end:
        return "representative_ORF"

    elif orf_end <= main_start:
        return "relative_upstream_ORF"

    elif orf_start >= main_end:
        return "relative_downstream_ORF"

    else:
        return "overlapping_ORF"

def analyze_transcript_orfs(transcript_id, db, bam_file, fasta_path):

    fasta_seq = load_transcript_sequence_from_genome(transcript_id, db, fasta_path)
    if fasta_seq is None:
        return []

    exon_regions, strand, seqid = get_transcript_features(db, transcript_id)
    if not exon_regions:
        return []

    orf_list = find_all_orf_candidates(
        fasta_seq, START_CODONS, STOP_CODONS, MIN_ORF_LENGTH_BP
    )

    if len(orf_list) == 0:
        return []

    # p-site periodicity calculation
    temp_results = []

    for orf_ts_start, orf_ts_end in orf_list:

        genomic_start = transcript_to_genomic_coordinate(exon_regions, orf_ts_start, strand)
        genomic_end = transcript_to_genomic_coordinate(exon_regions, orf_ts_end - 1, strand)

        if genomic_start is None or genomic_end is None:
            continue

        dna = fasta_seq[orf_ts_start:orf_ts_end]
        aa = str(Seq(dna).translate())[:-1]

        if aa.count('*') > 1:
            continue

        total_p_sites, ratio, frame0, frame_counts = calculate_periodicity_p_sites(orf_ts_start, orf_ts_end, bam_file, transcript_id, db)

        if (frame0 >= MIN_P_SITES_FOR_CODING) and (ratio >= MIN_FRAME_0_RATIO):

            temp_results.append({
                'transcript_id': transcript_id,
                'chromosome': seqid,
                'strand': strand,
                'genomic_start_loci': min(genomic_start, genomic_end),
                'genomic_end_loci': max(genomic_start, genomic_end),
                'orf_length_bp': orf_ts_end - orf_ts_start,
                'orf_status': 'Coding_Verified',
                'start_codon': dna[:3],
                'predicted_orf_p_sites': frame0,
                'total_orf_p_sites': total_p_sites,
                'frame_0_ratio (%)': ratio,
                'frame_0_p_sites': frame0,
                'frame_counts': frame_counts,
                'orf_dna_sequence': dna,
                'amino_acid_sequence': aa,
                'orf_ts_start': orf_ts_start,
                'orf_ts_end': orf_ts_end
            })

    if len(temp_results) == 0:
        return []

    # canonical flitration
    canonical_supported = [r for r in temp_results if r.get('orf_dna_sequence', '').startswith('ATG')]
    if len(canonical_supported) > 0:
        temp_results = canonical_supported

    main_orf = choose_main_orf(temp_results)

    main_start = main_orf["orf_ts_start"]
    main_end = main_orf["orf_ts_end"]

    final = []

    for r in temp_results:

        r["orf_class"] = classify_relative_orf(
            r["orf_ts_start"],
            r["orf_ts_end"],
            main_start,
            main_end
        )

        final.append(r)

    return final

def configure_psite_offsets(read_lengths_str, offsets_str):
    global P_SITE_OFFSETS

    lengths = [int(x) for x in read_lengths_str.split(',') if x.strip() != '']
    offsets = [int(x) for x in offsets_str.split(',') if x.strip() != '']

    if len(lengths) != len(offsets):
        raise ValueError('RPF read lengths and P-site offsets must have the same number of entries.')

    P_SITE_OFFSETS = dict(zip(lengths, offsets))
    print(f"Applied P-site offsets: {P_SITE_OFFSETS}")


def natural_chrom_key(chrom):
    chrom = str(chrom)
    c = chrom.replace('chr', '')

    if c.isdigit():
        return (0, int(c))
    if c == 'X':
        return (1, 23)
    if c == 'Y':
        return (1, 24)
    if c in ['M', 'MT']:
        return (1, 25)
    return (2, c)


# Main function
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python sORF_finder.py <sample_prefix> <RPF_Read_Lengths> <P_site_Offsets>')
        print('Example: python sORF_finder.py DLD1_alone 27,28,30 3,3,12')
        sys.exit(1)

    BAM_PREFIX = sys.argv[1]
    READ_LENGTHS_STR = sys.argv[2]
    OFFSETS_STR = sys.argv[3]
   
    try:
        configure_psite_offsets(READ_LENGTHS_STR, OFFSETS_STR)
    except ValueError as e:
        print(f"ERROR {e}")
        sys.exit(1)

    # input file paths
    DB_PATH = '/home/yewon/20_isoquant/tpm_percentage/GTF_majorisoform/GENCODEv34_addliftover_BIGv3_MajorIsoform.gtf.db'
    FASTA_PATH = '/home/yewon/Data_set/BRL/reference_genome/genome.fa'
    BAM_DIR = '/home/yewon/16_Ribo_seq/BRL_coculture/QCtoRiboTaper/7_STAR/MajorIsoforms/snoRNA_filtering/pcg_filtering/'
    BAM_PATH = f'{BAM_DIR}{BAM_PREFIX}.snoRNA_PCGfilt.Aligned.sortedByCoord.out.bam'
    GENE_LIST_PATH = '/home/yewon/16_Ribo_seq/BRL_coculture/RNAseq_readcounts/interest_genes.txt'
    TRANSCRIPT_TPM_PATH = "/home/yewon/20_isoquant/tpm_percentage/OUT.transcript_model/tpm5_or_10pct_majorisoform.tsv"

    print(f"BAM: {BAM_PATH}")
    print(f"Gene list: {GENE_LIST_PATH}")

    try:
        db = gffutils.FeatureDB(DB_PATH)
        bam_file = pysam.AlignmentFile(BAM_PATH, 'rb')
    except Exception as e:
        print(f"Failed to load DB or BAM: {e}")
        sys.exit(1)

    gene_list = []
    with open(GENE_LIST_PATH, 'r') as f:
        for line in f:
            item = line.strip()
            if item:
                gene_list.append(item)

    all_results = []

    for identifier in gene_list:
        print(f"\n--- Analyzing: {identifier} ---")

        if identifier.startswith('ENSG'):
            gene_ids_to_analyze = [identifier]
        else:
            gene_ids_to_analyze = get_gene_id_from_name(db, identifier)

        if not gene_ids_to_analyze:
            print(f"No gene_id found for: {identifier}")
            continue

        for gene_id in gene_ids_to_analyze:
            try:
                gene = db[gene_id]
                gene_name_to_use = get_gene_name_from_gene_feature(gene, identifier)
                transcript_ids = [t.id for t in db.children(gene, featuretype='transcript')]
            except gffutils.NoFeatError:
                print(f"gene_id not found in DB: {gene_id}")
                continue

            if not transcript_ids:
                print(f"No transcripts found for gene_id: {gene_id}")
                continue

            for tid in transcript_ids:
                transcript_results = analyze_transcript_orfs(tid, db, bam_file, FASTA_PATH)

                for result in transcript_results:
                    result['gene_id'] = gene_id
                    result['gene_name'] = gene_name_to_use
                    all_results.append(result)

            print(f"--- Done: {gene_id} ({gene_name_to_use}) ---")

    bam_file.close()

    df = pd.DataFrame(all_results)

    if df.empty:
        print('No ORFs passed filtering.')
        sys.exit(0)

    
    tpm_df = pd.read_csv(
    TRANSCRIPT_TPM_PATH,
    sep=r"\s+",
    header=None,
    names=["gene_name", "gene_id", "transcript_id", "TPM", "percentage"],
    engine="python")

    tpm_df["TPM"] = pd.to_numeric(tpm_df["TPM"], errors="coerce").fillna(0)
    tpm_df["percentage"] = pd.to_numeric(tpm_df["percentage"], errors="coerce").fillna(0)

    df = df.merge(
    tpm_df[["transcript_id", "TPM", "percentage"]],
    on="transcript_id",
    how="left")

    df["TPM"] = df["TPM"].fillna(0)
    df["percentage"] = df["percentage"].fillna(0)



    df['genomic_key'] = df.apply(
        lambda row: f"{row['chromosome']}:{row['genomic_start_loci']}-{row['genomic_end_loci']}:{row['strand']}",
        axis=1
    )

    df = df.sort_values(
    by=[
        'predicted_orf_p_sites',
        'frame_0_ratio (%)',
        'TPM',
        'percentage'
    ],
    ascending=[False, False, False, False])

    df = df.drop_duplicates(subset=['genomic_key'], keep='first').reset_index(drop=True)

    df['_chrom_sort'] = df['chromosome'].apply(natural_chrom_key)

    df = df.sort_values(
        by=['_chrom_sort', 'genomic_start_loci', 'genomic_end_loci', 'gene_name', 'transcript_id'],
        ascending=[True, True, True, True, True]
    ).reset_index(drop=True)


    column_order = [
    'chromosome',
    'gene_name',
    'gene_id',
    'transcript_id',
    'TPM',
    'percentage',
    'strand',
    'genomic_start_loci',
    'genomic_end_loci',
    'orf_length_bp',
    'orf_status',
    'start_codon',
    'predicted_orf_p_sites',
    'total_orf_p_sites',
    'frame_0_ratio (%)',
    'orf_class',
    'frame_counts',
    'orf_dna_sequence',
    'amino_acid_sequence']

    for col in column_order:
        if col not in df.columns:
            df[col] = ''

    df_out = df[column_order]

    output_path = f'{BAM_PREFIX}_sORF_results.tsv'
    df_out.to_csv(output_path, sep='\t', index=False)

    print(f"\nSaved TSV output: {output_path}")
    print(f"Number of reported ORFs: {len(df_out)}")
