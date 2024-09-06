from Bio  import SeqIO
def fq_from_abi(file_path: str) -> None:
    '''读取ab1文件并将序列信息存入fq格式字符串'''
    fq_str = ""
    for seq in SeqIO.parse(file_path, "abi"):
        seq_id = seq.name
        # base signal
        abif_raw = seq.annotations["abif_raw"]
        data_a = list(abif_raw["DATA9"])
        data_c = list(abif_raw["DATA10"])
        data_g = list(abif_raw["DATA11"])
        data_t = list(abif_raw["DATA12"])
        for a,g,c,t in zip(data_a,data_g,data_c,data_t):
            print(f"{a}\t{t}\t{c}\t{g}")
        # signal qc

        # trim seq and qual,storage in dict
        # trimmed_seq, trimmed_qual = trim_static(seq)
        # trimmed qual control
        # if not trimmed_qual_qc(trimmed_qual):
        #     print(f"{file_path} qc unpass !!!")
        #     return ''
        # trimmed_qual_str = "".join(trimmed_qual)
        # fq_str = f"@{seq_id}\n{trimmed_seq}\n+\n{trimmed_qual_str}"
    return fq_str


fq_from_abi("B103-(T3389)pUp-pDown-flank-R.ab1")