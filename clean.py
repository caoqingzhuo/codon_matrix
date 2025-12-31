from extract_cds import build_clean_cds

build_clean_cds(
    input_path="./p_pacificus.PRJNA12644.WS285.CDS_transcripts.fa",
    input_format="fasta",
    out_fasta="./p_pacificus_clean_CDS.fasta",
    table_id=1,
    min_len=90
)
