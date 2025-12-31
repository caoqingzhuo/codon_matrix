from Bio import SeqIO
from Bio.Data import CodonTable


def get_start_stop_codons(table_id):
    table = CodonTable.unambiguous_dna_by_id[table_id]
    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)
    return start_codons, stop_codons


def clean(seq, table_id=1, min_len=90):
    s = str(seq).upper().replace("U", "T")
    start_codons, stop_codons = get_start_stop_codons(table_id)

    start_pos = None
    for i in range(0, len(s) - 2):
        if s[i:i+3] in start_codons:
            start_pos = i
            break
    if start_pos is None:
        return None

    for i in range(start_pos, len(s) - 2, 3):
        if s[i:i+3] in stop_codons:
            cds = s[start_pos:i+3] 
            if len(cds) % 3 != 0:
                return None
            if len(cds) < min_len:
                return None
            return cds
    return None


def extract_sequences_from_genbank(gb_path):
    record = SeqIO.read(gb_path, "genbank")
    out = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        cds_seq = feat.extract(record.seq)
        gene = feat.qualifiers.get("gene", ["unknown_gene"])[0]
        product = feat.qualifiers.get("product", ["unknown_product"])[0]
        header = "%s|gene=%s|product=%s" % (record.id, gene, product)
        out.append((header, cds_seq))
    return out


def extract_sequences_from_fasta(fasta_path):
    out = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        out.append((rec.id, rec.seq))
    return out


def write_fasta(records, out_fasta):
    with open(out_fasta, "w") as f:
        for header, seq in records:
            f.write(">%s\n" % header)
            s = str(seq)
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + "\n")


def build_clean_cds(input_path, input_format, out_fasta, table_id=1, min_len=90):
    if input_format == "genbank":
        raw_records = extract_sequences_from_genbank(input_path)
    elif input_format == "fasta":
        raw_records = extract_sequences_from_fasta(input_path)

    clean_records = []
    for header, seq in raw_records:
        cds = clean(seq, table_id=table_id, min_len=min_len)
        if cds is None:
            continue
        clean_records.append((header + "|cleaned_len=%d" % len(cds), cds))
    write_fasta(clean_records, out_fasta)
