from Bio import SeqIO

def peak_evaluation():
    ...


def trim_seq_by_qual(seq_record, threshold=20):
    # 获取序列的质量分数
    qual_scores = seq_record.letter_annotations["phred_quality"]
    
    # 找到第一个质量低于阈值的位置
    trim_position = len(qual_scores)
    for i, score in enumerate(qual_scores):
        if score < threshold:
            trim_position = i
            break
    
    # 修剪序列和质量分数
    seq_record.seq = seq_record.seq[:trim_position]
    seq_record.letter_annotations["phred_quality"] = qual_scores[:trim_position]
    return seq_record


def parse_abi_file(file_path:str) -> None:
    for seq in SeqIO.parse(file_path,"abi"):
        seq_id, sequence, qual_array = seq.name, seq.seq, seq.letter_annotations['phred_quality']
        abif_raw = seq.annotations["abif_raw"]
        data_g = list(abif_raw["DATA9"])
        data_a = list(abif_raw["DATA10"])
        data_t = list(abif_raw["DATA11"])
        data_c = list(abif_raw["DATA12"])
        trimmed_record = trim_seq_by_qual(seq)
        print(trimmed_record.seq)
        print(trimmed_record.letter_annotations['phred_quality'])


def main() -> None:
    input_file = "/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1"
    parse_abi_file(input_file)
    

if __name__ == "__main__":
    main()