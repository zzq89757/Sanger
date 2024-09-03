from collections import defaultdict
from Bio import SeqIO
import subprocess
import pandas as pd
from pysam import AlignedSegment



    
def recognized_well_by_file_name(file_name:str = "HA-1-4-B-6") -> int:
    '''根据文件名对应的子板及孔位获取拆分前的孔号'''

    sub_dict = defaultdict(list)
    
    for i in range(1, 17):
        for j in range(1, 25):
            if i % 2 == 1:
                if j % 2 == 1:
                    sub_dict['1'].append((i - 1) * 24 + j)
                else:
                    sub_dict['2'].append((i - 1) * 24 + j)
            else:
                if j % 2 == 1:
                    sub_dict['3'].append((i - 1) * 24 + j)
                else:
                    sub_dict['4'].append((i - 1) * 24 + j)
    

    subplate, raw_num, col_num = file_name.split("-")[-3:]
    raw_num = ord(raw_num) - 65
    sub_well_idx = raw_num * 12 + int(col_num) - 1
    well = sub_dict[subplate][sub_well_idx]
    return well
    
        


def generate_ref(sg_table:str, vector_seq:str, ref_path:str) -> defaultdict:
    '''根据sgRNA序列表生成每个plate->well对应的reference'''
    # read sg seq table
    sg_df = pd.read_excel(sg_table)

    # extracat raw vector seq
    vector_seq:str
    for v in SeqIO.parse(vector_seq,'fasta'):vector_seq = str(v.seq)
    vector_li = vector_seq.split("N" * 20)
    
    # generate reference by sg df and raw vector seq
    idx = 1
    tmp_li = [''] * 5
    for plate_str, well, sgn_seq, ori in zip(sg_df['Plate #'], sg_df['Well #'], sg_df['sgRNA sequence'], sg_df['Strand']):
        plate, sg_idx = plate_str.split("_sg")
        tmp_li[int(sg_idx) - 1] = vector_li[int(sg_idx) - 1] + sgn_seq
        if idx % 4 == 0:
            tmp_li[4] = vector_li[4]
            # storage finished,join and output as reference file
            well_ref = ''.join(tmp_li)
            out_ref = open(f"{ref_path}/{plate}_{well}.fa",'w')
            out_ref.write(f">{plate}_{well}\n{well_ref}\n")
            # construct index
            subprocess.run(f"bwa index {ref_path}/{plate}_{well}.fa", shell=True)
            tmp_li = [''] * 5
        idx += 1  
    
    

    

def peak_evaluation():
    '''检测下机数据重叠峰'''
    ...

def trim_static(seq_record, start:int=50, end:int=800):
    '''裁剪下机数据，固定保留50-800部分'''
    # return seq_record.seq[50:800], seq_record.letter_annotations["phred_quality"][50:800]
    return seq_record.seq[50:800], "F"*750
    


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


def fq_from_abi(file_path:str) -> defaultdict:
    '''读取ab1文件并将序列信息存入fq格式字符串'''
    fq_str = ""
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
        fq_str = f"@{seq_id}\n{trimmed_seq}\n+\n{trimmed_qual}"
    return fq_str



def parse_alignment_result(aln_res_str:str):
    a = AlignedSegment()
    a.flag = 4
    a.cigarstring = "47M150S"
    print(a.get_cigar_stats())
    
    

def process_align(ref_file:str, fq_str:str) -> None:
    '''将下机数据比对到参考并处理结果'''
    # alignment by bwa and fetch alignment result
    res = subprocess.Popen(f"echo -e \"{fq_str}\" | bwa mem -t 24 {ref_file} /dev/stdin",shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    aln_res_str = str(res.stdout.read(),"utf-8").rstrip().split("\n")[-1]
    # process alignment result
    parse_alignment_result(aln_res_str)
    

    

def main() -> None:
    input_file = "/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1"
    # generate_ref("/home/wayne/Project/SC/Sanger/Ho_tf1.xlsx","/home/wayne/Project/SC/Sanger/raw_vector.fa","./newref/")
    well = recognized_well_by_file_name()
    ref_file = f"{plate}_{well}.fa"
    fq_str = fq_from_abi(input_file)
    process_align("raw_vector.fa", fq_str)
    

if __name__ == "__main__":
    main()