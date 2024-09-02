from collections import defaultdict
from Bio import SeqIO,Align


def parse_fq_ref_map_table() -> defaultdict:
    '''将下机数据与其对应参考的映射关系存入字典'''

def generate_ref(vector_seq:str, sg_table:str):
    '''根据原始载体和sgRNA序列生成不同的参考文件'''
    ...

def peak_evaluation():
    '''检测下机数据重叠峰'''
    ...

def trim_static(seq_record, start:int=50, end:int=800):
    '''裁剪下机数据，固定保留50-800部分'''
    return seq_record.seq[50:800],seq_record.letter_annotations["phred_quality"][50:800]
    


def trim_seq_by_qual(seq_record, threshold=20):
    '''根据质量值截取序列'''
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


def data_from_abi(file_path:str) -> defaultdict:
    '''读取ab1文件并将序列信息存入字典'''
    data_dict = defaultdict(lambda:defaultdict())
    for seq in SeqIO.parse(file_path,"abi"):
        seq_id= seq.name
        # base signal
        abif_raw = seq.annotations["abif_raw"]
        data_g = list(abif_raw["DATA9"])
        data_a = list(abif_raw["DATA10"])
        data_t = list(abif_raw["DATA11"])
        data_c = list(abif_raw["DATA12"])
        # trimmed_record = trim_seq_by_qual(seq)
        # trim seq and qual,storage in dict
        trimmed_seq, trimmed_qual = trim_static(seq)
        data_dict[seq_id]['seq'] = trimmed_seq
        data_dict[seq_id]['qual'] = trimmed_qual
    return data_dict

def align2ref(ref_file:str, data_dict:defaultdict) -> None:
    '''构建参考序列索引并将下机数据比对到参考'''
    from os import system
    system(f"bwa index {ref_file}")
    ...



def main() -> None:
    input_file = "/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1"
    data_dict = data_from_abi(input_file)
    print(data_dict)
    

if __name__ == "__main__":
    main()