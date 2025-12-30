from Bio import SeqIO
from Bio.Seq import Seq

gb_path = "./sequence.gb"          
out_fasta = "mito_CDS.fasta"      

record = SeqIO.read(gb_path, "genbank")

print("Record ID:", record.id)
print("Length:", len(record.seq))
print("Description:", record.description)
print("Number of features:", len(record.features))

cds_seqs = []
for feat in record.features:
    if feat.type != "CDS":
        continue

    # 取出CDS序列（Biopython会自动处理正负链）
    cds: Seq = feat.extract(record.seq)

    # 一些注释信息（尽量取到一个能标识的名字）
    gene = feat.qualifiers.get("gene", ["unknown_gene"])[0]
    product = feat.qualifiers.get("product", ["unknown_product"])[0]

    s = str(cds).upper().replace("U", "T")
    if any(b not in "ATCG" for b in s):
        continue
    if len(s) < 30:
        continue    

    header = f"{record.id}|gene={gene}|product={product}|len={len(s)}"
    cds_seqs.append((header, s))

print("Extracted CDS count:", len(cds_seqs))

# 写出fasta
with open(out_fasta, "w") as f:
    for h, seq in cds_seqs:
        f.write(f">{h}\n")
        # 每行60字符
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

print("Wrote:", out_fasta)
