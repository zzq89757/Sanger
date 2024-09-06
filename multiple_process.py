from collections import defaultdict, deque
from math import ceil
from pathlib import Path
import subprocess
from Bio import SeqIO
from pysam import AlignmentFile
import pandas as pd


def cut_seq_obtain_pos(vector_seq: str) -> list:
    '''根据原始载体构建sgRNA位置列表以及载体不含sgRNA的切片列表'''
    for v in SeqIO.parse(vector_seq, 'fasta'):
        vector_seq = str(v.seq)
    vector_li = vector_seq.split("N" * 20)
    # record sgRNA pos by strlen
    sg_pos_li = []
    front_length = 0
    for i, v in enumerate(vector_li):
        if i == 4:
            continue
        pos = len(v)
        front_length += pos
        sg_pos_li.append(front_length)
        front_length += 20

    return sg_pos_li, vector_li


def generate_ref(sg_table: str, vector_li: list, ref_path: str) -> None:
    '''根据sgRNA序列表生成每个plate->well对应的reference'''
    # read sg seq table
    sg_df = pd.read_excel(sg_table)

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
            out_ref = open(f"{ref_path}/{plate}_{well}.fa", 'w')
            out_ref.write(f">{plate}_{well}\n{well_ref}\n")
            # construct index
            subprocess.run(
                f"bwa index {ref_path}/{plate}_{well}.fa", shell=True)
            tmp_li = [''] * 5
        idx += 1


# 提取同一个well中的多个ab1文件中的序列信息 保存为fq后 bwa比对 处理bam文件


def recognized_well_by_file_name(file_name: str = "HA-1-4-B-6") -> int:
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


def classify_file_by_well(file_path: Path) -> defaultdict:
    '''输入下机数据路径,将其分配至对应的well'''
    well_fq_file_dict = defaultdict(list)
    file_li = file_path.glob("*.ab*")
    for file in file_li:
        file_name = file.split("-")[:5]
        well = recognized_well_by_file_name(file_name)
        well_fq_file_dict[well].append(file)
    return well_fq_file_dict


def set_end_pos(start_pos: int, next_base_pos: int, end_pos: int, trace_len: int) -> int:
            if next_base_pos:
                end_pos = start_pos + ceil((next_base_pos - start_pos) / 2)
            else:
                end_pos = trace_len
            return end_pos
        


def peak_qc(a_trace, g_trace, t_trace, c_trace) -> bool:
    '''对每个轨道进行信号重叠以及峰型分析'''
    
    

def trim_static(seq_record, start: int = 50, end: int = 800) -> list:
    '''裁剪下机数据，固定保留50-800部分'''
    # return seq_record.seq[50:800], "F" * 750
    return [seq_record.seq[50:800], seq_record.letter_annotations["phred_quality"][50:800]]


def trimmed_qual_qc(qc_array) -> bool:
    # calc mean qual

    # calc q20 and q30
    q_20 = len([x for x in qc_array if x >= 20])/len(qc_array)
    q_30 = len([x for x in qc_array if x >= 30])/len(qc_array)
    if q_20 < 0.9 or q_30 < 0.8:
        return False
    return True


def fq_from_abi(file_path: str) -> defaultdict:
    '''读取ab1文件并将序列信息存入fq格式字符串'''
    fq_str = ""
    for seq in SeqIO.parse(file_path, "abi"):
        seq_id = seq.name
        # base signal
        abif_raw = seq.annotations["abif_raw"]
        g_trace = list(abif_raw["DATA9"])
        a_trace = list(abif_raw["DATA10"])
        t_trace = list(abif_raw["DATA11"])
        c_trace = list(abif_raw["DATA12"])
        # signal qc

        # trim seq and qual,storage in dict
        trimmed_seq, trimmed_qual = trim_static(seq)
        # trimmed qual control
        if not trimmed_qual_qc(trimmed_qual):
            print(f"{file_path} qc unpass !!!")
            return ''
        trimmed_qual_str = "".join(trimmed_qual)
        fq_str = f"@{seq_id}\n{trimmed_seq}\n+\n{trimmed_qual_str}"
    return fq_str


def extract_data(well_fq_file_dict: defaultdict, output_fq_path: str) -> None:
    '''提取每个well中的序列信息并将属于同一个well的序列存为fq'''
    for well, file_li in well_fq_file_dict.items():
        output_handle = open(f"{output_fq_path}/{well}_trimed.fq", 'w')
        for file in file_li:
            fq_str = fq_from_abi(file)
            output_handle.write(fq_str)


def sgRNA_detective(start: int, end: int, sg_pos_li: list) -> list[int]:
    '''比对结果覆盖到了第几条sgRNA'''
    sgRNA_cover_idx_li = []
    for idx, pos in enumerate(sg_pos_li):
        if end > pos + 20 and pos > start:
            sgRNA_cover_idx_li.append(idx)

    return sgRNA_cover_idx_li

def coverage_check(cov_start: int, cov_end: int, aln_region_li:list) -> bool:
    '''检查是否完全覆盖目标区域'''
    cov_dict = defaultdict(0)
    for reg in aln_region_li:
        for i in range(reg[0] + 1, reg[1] + 1):
            cov_dict[i] = 1
    for i in range(cov_start + 1, cov_end + 1):
        if not cov_dict[i]:return False
    return True
    

def parse_alignment_result(output_bam: str, sg_pos_li: list, well_qc_dict: defaultdict) -> None:
    '''处理比对的结果'''
    well = int(output_bam.split("-")[0])
    well_qc_dict[well]['sgRNA'] = [0] * 4
    well_qc_dict[well]['coverage'] = 0
    aln_region_li = []
    for aln in AlignmentFile(output_bam, 'r', threads=16):
        # if with mismath or softclip warnning
        if sum(aln.get_cigar_stats()[0][1:]):
            print(f"{aln.query_name} has mismatch or softclip !!!")
            well_qc_dict[well]['invalid_seq_num'].append(aln.query_name)
            return 0
        # sgRNA detective
        cover_idx_li: list[int] = sgRNA_detective(
            aln.reference_start, aln.reference_end, sg_pos_li)
        for i in cover_idx_li:
            well_qc_dict[well]['sgRNA'][i] = 1
        aln_region_li.append([aln.reference_start, aln.reference_end])
    if coverage_check(1, 1000, aln_region_li):well_qc_dict[well]['coverage'] = 1


def process_alignment(ref_file: str, input_fq: str, output_bam: str, sg_pos_li: list, well_qc_dict: defaultdict) -> None:
    '''执行bwa命令后处理bam文件'''
    # alignment by bwa and fetch alignment result
    res = subprocess.Popen(f"bwa mem -t 24 {ref_file} {input_fq} > {output_bam}",
                           shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # process alignment result
    parse_alignment_result(output_bam, sg_pos_li, well_qc_dict)


def qc_dict_to_table(well_qc_dict: defaultdict, output_table: str) -> None:
    '''将包含qc信息的dict转为表格并存储为文件'''
    out_table_handle = open(output_table, 'w')
    out_table_handle.write("Well\tsgRNA1\tsgRNA2\tsgRNA3\tsgRNA4\tall\tcoverage\n")
    for well in well_qc_dict.keys():
        all_detective = not 0 in well_qc_dict[well]['sgRNA']
        well_qc_str = str(well) + "\t" + "\t".join(
            [str(x) for x in well_qc_dict[well]['sgRNA']]) + "\t" + str(all_detective) + "\t" + str(well_qc_dict[well]['coverage']) + "\n"
        out_table_handle.write(well_qc_str)


def main() -> None:
    # construct reference of each well
    sg_pos_li, vector_li = cut_seq_obtain_pos(
        "/home/wayne/Project/SC/Sanger/Ho_tf1.xlsx", "/home/wayne/Project/SC/Sanger/raw_vector.fa")
    generate_ref("/home/wayne/Project/SC/Sanger/Ho_tf1.xlsx",
                 vector_li, "./newref/")
    # sanger file process
    file_dict = classify_file_by_well("/home/wayne/Project/SC/Sanger/")
    output_fq_path = "./fq/"
    extract_data(file_dict, output_fq_path)
    well_li = [str(x).split("_")[0] for x in Path(output_fq_path).glob("*fq")]
    well_qc_dict = defaultdict(lambda: defaultdict(list))
    # alignment and output
    for well in well_li:
        ref_file = ""
        input_fq = ""
        output_bam = ""
        process_alignment()
    qc_dict_to_table(well_qc_dict, "./res.tsv")


if __name__ == "__main__":
    main()
