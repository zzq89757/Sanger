from collections import defaultdict, deque
from pathlib import Path
import subprocess
from Bio import SeqIO
from pysam import AlignmentFile
import pandas as pd



def cut_seq_obtain_pos(vector_seq:str) -> list:
    '''根据原始载体构建sgRNA位置列表以及载体不含sgRNA的切片列表'''
    for v in SeqIO.parse(vector_seq,'fasta'):vector_seq = str(v.seq)
    vector_li = vector_seq.split("N" * 20)
    # record sgRNA pos by strlen  
    sg_pos_li = []
    front_length = 0
    for i,v in enumerate(vector_li):
        if i == 4:continue
        pos = len(v) 
        front_length += pos
        sg_pos_li.append(front_length)
        front_length += 20
        
        
    return sg_pos_li , vector_li


def generate_ref(sg_table:str, vector_li:list, ref_path:str) -> None:
    '''根据sgRNA序列表生成每个plate->well对应的reference'''
    # read sg seq table
    sg_df = pd.read_excel(sg_table)

    
    # sg_pos_li , vector_li = cut_seq_obtain_pos(vector_seq)
    
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
            # out_ref = open(f"{ref_path}/{plate}_{well}.fa",'w')
            # out_ref.write(f">{plate}_{well}\n{well_ref}\n")
            # construct index
            # subprocess.run(f"bwa index {ref_path}/{plate}_{well}.fa", shell=True)
            tmp_li = [''] * 5
        idx += 1 
    



# 提取同一个well中的多个ab1文件中的序列信息 保存为fq后 bwa比对 处理bam文件

    
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
   


def classify_file_by_well(file_path:Path) -> defaultdict:
    '''输入下机数据路径,将其分配至对应的well'''
    well_fq_file_dict = defaultdict(list)
    file_li = file_path.glob("*.ab*")
    for file in file_li:
        file_name = file.split("-")[:5]
        well = recognized_well_by_file_name(file_name)
        well_fq_file_dict[well].append(file)
    return well_fq_file_dict




def trim_static(seq_record, start:int=50, end:int=800):
    '''裁剪下机数据，固定保留50-800部分'''
    # return seq_record.seq[50:800], seq_record.letter_annotations["phred_quality"][50:800]
    return seq_record.seq[50:800], "F" * 750
    


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


        



def extract_data(well_fq_file_dict:defaultdict) -> None:
    '''提取每个well中的序列信息并将属于同一个well的序列存为fq'''
    for well, file_li in well_fq_file_dict.items():
        output_handle = open(f"{well}_trimed.fq",'w')
        for file in file_li:
            fq_str = fq_from_abi(file)
            output_handle.write(fq_str)
    
    
def detective_alignment_result():
    '''检测是否存在异常的结果'''
    ...
    
    
def process_alignment(ref_file:str, input_fq:str, output_bam:str, sg_pos_li):
    '''执行bwa命令后处理bam文件'''
    # alignment by bwa and fetch alignment result
    res = subprocess.Popen(f"bwa mem -t 24 {ref_file} {input_fq} > ",shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    aln_res_str = str(res.stdout.read(),"utf-8").rstrip().split("\n")[-1]
    sg_pos_li = generate_ref("/home/wayne/Project/SC/Sanger/Ho_tf1.xlsx","/home/wayne/Project/SC/Sanger/raw_vector.fa","./newref/")
    # process alignment result
    detective_alignment_result(aln_res_str, sg_pos_li)



def main():
    # construct reference of each well
    sg_pos_li, vector_li = cut_seq_obtain_pos("/home/wayne/Project/SC/Sanger/Ho_tf1.xlsx","/home/wayne/Project/SC/Sanger/raw_vector.fa")
    generate_ref("/home/wayne/Project/SC/Sanger/Ho_tf1.xlsx",vector_li,"./newref/")
    # process alignment
    file_dict = classify_file_by_well("/home/wayne/Project/SC/Sanger/")
    input_file = "/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1"
    well = recognized_well_by_file_name()
    plate = 0
    ref_file = f"{plate}_{well}.fa"
    fq_str = fq_from_abi(input_file)
    process_align("raw_vector.fa", fq_str, sg_pos_li)