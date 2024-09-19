from collections import defaultdict
from os import system
from pathlib import Path
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment
import pandas as pd


def cut_seq_obtain_pos(seq_path: str) -> list:
    """根据原始载体构建sgRNA位置列表以及载体不含sgRNA的切片列表和元件位置字典"""
    vector_seq: str = ""
    feature_dict = {}
    for record in SeqIO.parse(seq_path, "genbank"):
        for feature in record.features:
            start = int(feature.location.start)
            end = int(feature.location.end)

            # 如果 feature 有 /label，则存入字典
            if "label" in feature.qualifiers:
                label = feature.qualifiers["label"][0]
                # 将位置范围作为键，label 作为值存入字典
                feature_dict[(start, end)] = label
        vector_seq = str(record.seq)

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

    return sg_pos_li, vector_li, feature_dict


def find_label_by_position(feature_dict: dict, position: int) -> str:
    """
    根据输入的位置从字典中查找对应的 /label。

    :param feature_dict: feature 位置信息和 /label 的字典
    :param position: 输入的位置（从 0 开始）
    :return: 对应位置的 label，若未找到则返回 "No label found."
    """
    for (start, end), label in feature_dict.items():
        if start <= position <= end:
            return label
    return "-"


def generate_ref(sg_table: str, vector_li: list, ref_path: str) -> defaultdict:
    """根据sgRNA序列表生成每个plate->well对应的reference"""
    well_ref_dict = defaultdict(str)
    # aln sg seq table
    sg_df = pd.read_excel(sg_table)

    # generate reference by sg df and raw vector seq
    idx = 1
    tmp_li = [""] * 5
    for plate_str, well, sgn_seq, ori in zip(
        sg_df["Plate #"], sg_df["Well #"], sg_df["sgRNA sequence"], sg_df["Strand"]
    ):
        plate, sg_idx = plate_str.split("_sg")
        plate: str = plate.upper().replace('TF','').replace('_','-')
        tmp_li[int(sg_idx) - 1] = vector_li[int(sg_idx) - 1] + sgn_seq
        if idx % 4 == 0:
            tmp_li[4] = vector_li[4]
            # storage finished,join and output as reference file
            well_ref = "".join(tmp_li)
            well_ref_dict[well] = well_ref
            if not Path(f"{ref_path}/{plate}_{well}.fa").exists():
                # write into fastq file
                out_ref = open(f"{ref_path}/{plate}_{well}.fa", "w")
                out_ref.write(f">{plate}_{well}\n{well_ref}\n")
                # construct index
                system(f"bwa index {ref_path}/{plate}_{well}.fa")
            tmp_li = [""] * 5
        idx += 1
    return well_ref_dict





def recognized_well_by_file_name(file_name: str = "HA-1-4-A01-P1_A01#2.ab1") -> int:
    """根据文件名对应的子板及孔位获取拆分前的孔号"""

    sub_dict = defaultdict(list)
    for i in range(1, 17):
        for j in range(1, 25):
            if i % 2 == 1:
                if j % 2 == 1:
                    sub_dict["1"].append((i - 1) * 24 + j)
                else:
                    sub_dict["2"].append((i - 1) * 24 + j)
            else:
                if j % 2 == 1:
                    sub_dict["3"].append((i - 1) * 24 + j)
                else:
                    sub_dict["4"].append((i - 1) * 24 + j)

    # subplate, raw_num, col_num = file_name.split("-")[-3:]
    subplate, raw_col = file_name.split("#")[0].split("_")[0].split("-")[2:4]
    raw_num = raw_col[0]
    col_num = int(raw_col[1:])
    raw_num = ord(raw_num) - 65
    sub_well_idx = raw_num * 12 + int(col_num) - 1
    well = sub_dict[subplate][sub_well_idx]
    return well


def trim_static(seq_record, start: int = 50, end: int = 800) -> list:
    """裁剪下机数据,固定保留50-800部分"""
    # return seq_record.seq[50:800], "F" * 750
    return [
        seq_record.seq[start:end],
        seq_record.letter_annotations["phred_quality"][start:end],
    ]


def trimmed_qual_qc(qc_array) -> bool:
    # calc mean qual

    # calc q20 and q30

    q_20 = len([x for x in qc_array if x >= 20]) / len(qc_array)
    q_30 = len([x for x in qc_array if x >= 30]) / len(qc_array)
    if q_20 < 0.9 or q_30 < 0.8:
        return False
    return True


def fq_from_abi(
    trim_start: int, trim_end: int, well_qc_dict: defaultdict, file_path: str
) -> defaultdict:
    """读取ab1文件并将序列信息存入fq格式字符串"""
    fq_str = ""
    for seq in SeqIO.parse(file_path, "abi"):
        prefix_li = str(seq.name).split("_")[0].split("-")
        subplate = "-".join(prefix_li[:4])
        pcr_cycle = '' if len(prefix_li) < 6 else f"#{prefix_li[5]}"
        sub_cycle_info = subplate + pcr_cycle
        seq_id = seq.name
        # base signal
        abif_raw = seq.annotations["abif_raw"]
        g_trace = list(abif_raw["DATA9"])
        a_trace = list(abif_raw["DATA10"])
        t_trace = list(abif_raw["DATA11"])
        c_trace = list(abif_raw["DATA12"])
        # signal qc

        # trim seq and qual,storage in dict
        trimmed_seq, trimmed_qual = trim_static(seq, trim_start, trim_end)
        if not trimmed_qual:
            well_qc_dict[sub_cycle_info]["qc_failed"] = [1]
            print(f"{sub_cycle_info} no signal !!!")
            return ""
        # trimmed qual control
        if not trimmed_qual_qc(trimmed_qual):
            well_qc_dict[sub_cycle_info]["qc_failed"] = [1]
            print(f"{sub_cycle_info} qc unpass !!!")
            return ""
        trimmed_qual_str = "".join(["F" for x in trimmed_qual])
        fq_str = f"@{seq_id}\n{trimmed_seq}\n+\n{trimmed_qual_str}\n"
    return fq_str


def generate_fq(
    trim_start: int,
    trim_end: int,
    well_qc_dict: defaultdict,
    input_path: str,
    output_fq_path: str,
) -> None:
    """提取每个well中的序列信息并将属于同一个well的序列存为fq"""
    file_li = Path(input_path).glob("*.ab*")
    for file in file_li:
        prefix_li = file.name.split("_")[0].split("-")
        subplate = "-".join(prefix_li[:4])
        pcr_cycle = '' if len(prefix_li) < 6 else f"#{prefix_li[5]}"
        sub_cycle_info = subplate + pcr_cycle
        well = recognized_well_by_file_name(subplate)
        output_handle = open(f"{output_fq_path}/{sub_cycle_info}_{well}.fq", "a")
        fq_str = fq_from_abi(trim_start, trim_end, well_qc_dict, file)
        output_handle.write(fq_str)


def stop_codon_check(sub_cycle_info: str, segment: str, well_qc_dict: defaultdict) -> None:
    """检测突变后的碱基是否导致终止密码子的产生"""
    stop_codon_li = ["UAA", "UAG", "UGA"]
    # stop codon detective by ORF start and SNP position

    # if not in cds, continue

    # in cds, get ORF and codon after mutation
    for i in range(3):
        codon = segment[i : i + 3]
        # in stop codon li and in cds
        if codon in stop_codon_li:
            well_qc_dict[sub_cycle_info]["stop_codon_from_snp"] = [1]


def mismatch_check(
    sub_cycle_info: str,
    aln: AlignedSegment,
    well_qc_dict: defaultdict,
    detective_end: int,
    feature_dict: dict,
) -> bool:
    """检查错配数目,等于1报点,大于1报错"""
    all_snp_count = sum(aln.get_cigar_stats()[0][1:])
    if all_snp_count > 1:
        print(f"{aln.query_name} has mismatch or softclip !!!")
        well_qc_dict[sub_cycle_info]["indel_soft"] += [1]
        return True
    elif all_snp_count == 1:
        # get snp position
        md_tag = aln.get_tag("MD")
        md_snp_idx = (
            md_tag.find("A") + 1
            or md_tag.find("G") + 1
            or md_tag.find("C") + 1
            or md_tag.find("T") + 1
        )
        forward_len = int(md_tag[: md_snp_idx - 1])
        pos = aln.reference_start + forward_len + 1
        component = find_label_by_position(feature_dict, pos)
        mismatch_str = f"{pos}<{component}>:{md_tag[md_snp_idx - 1]}->{aln.query_alignment_sequence[forward_len]}"
        if pos < detective_end and mismatch_str not in well_qc_dict[sub_cycle_info]["mismatch"]:
            # stop codon check
            # print(aln.query_name)
            stop_codon_check(
                sub_cycle_info,
                aln.query_alignment_sequence[forward_len - 2 : forward_len + 3],
                well_qc_dict,
            )
            well_qc_dict[sub_cycle_info]["mismatch"].append(
                f"{pos}<{component}>:{md_tag[md_snp_idx -1 ]}->{aln.query_alignment_sequence[forward_len]}"
            )
    # else:
    #     well_qc_dict[well]['mismatch'] = '0'
    return False


def sgRNA_detective(
    aln: AlignedSegment, sg_pos_li: list
) -> list[int]:
    """比对结果覆盖到了第几条sgRNA"""
    sgRNA_cover_idx_li = []
    q_start = aln.reference_start
    for idx, pos in enumerate(sg_pos_li):
        if (
            aln.query_sequence[pos - q_start : pos - q_start + 20].upper()
            == aln.get_reference_sequence()[pos : pos + 20].upper()
        ):
            sgRNA_cover_idx_li.append(idx)

    return sgRNA_cover_idx_li


def coverage_check(cov_start: int, cov_end: int, aln_region_li: list) -> bool:
    """检查是否完全覆盖目标区域"""
    cov_dict = defaultdict(int)
    for reg in aln_region_li:
        for i in range(reg[0] + 1, reg[1] + 1):
            cov_dict[i] = 1
    for i in range(cov_start + 50, cov_end - 50):
        if not cov_dict[i]:
            return False
    return True


def parse_alignment_result(
    output_bam: str,
    sg_pos_li: list,
    well_qc_dict: defaultdict,
    feature_dict: dict,
) -> None:
    """处理比对的结果"""
    prefix_li = Path(output_bam).name.split(".")[0].split("_")
    sub_cycle_info, well = prefix_li
    well_qc_dict[sub_cycle_info]["sgRNA"] = [0] * 4
    well_qc_dict[sub_cycle_info]["coverage"] = 0
    aln_region_li = []
    for aln in AlignmentFile(output_bam, "r", threads=16):
        # skip supplementary and secondary
        if aln.is_supplementary or aln.is_secondary:
            continue
        # mismatch check
        if mismatch_check(sub_cycle_info, aln, well_qc_dict, sg_pos_li[-1] + 100, feature_dict):
            continue
        # sgRNA detective
        cover_idx_li: list[int] = sgRNA_detective(aln, sg_pos_li)
        for i in cover_idx_li:
            well_qc_dict[sub_cycle_info]["sgRNA"][i] = 1
        aln_region_li.append([aln.reference_start, aln.reference_end])
    # coverage check after all fastq mapped to ref
    if coverage_check(5914, sg_pos_li[-1] + 150, aln_region_li):
        well_qc_dict[sub_cycle_info]["coverage"] = 1


def process_alignment(
    ref_file: str,
    input_fq: str,
    output_bam: str,
    sg_pos_li: list,
    well_qc_dict: defaultdict,
    feature_dict: dict,
) -> None:
    """执行bwa命令后处理bam文件"""
    # alignment by bwa and fetch alignment result
    if not Path(output_bam).exists():
        # print(f"bwa mem -t 24 {ref_file} {input_fq} > {output_bam}")
        system(f"bwa mem -t 24 {ref_file} {input_fq} > {output_bam}")

    # process alignment result
    parse_alignment_result(
        output_bam, sg_pos_li, well_qc_dict, feature_dict
    )


def qc_dict_to_table(
    well_qc_dict: defaultdict, output_table: str
) -> None:
    """将包含qc信息的dict转为表格并存储为文件"""
    out_table_handle = open(output_table, "w")
    out_table_handle.write(
        "Subplate_well\tWell\tsgRNA1\tsgRNA2\tsgRNA3\tsgRNA4\tall_sgRNA\tcoverage\tqc_failed\tmismatch\tindel_soft\n"
    )
    # print(well_qc_dict.keys())
    # for i in well_qc_dict.keys():
    #     print(i)
    #     recognized_well_by_file_name(i)
    # exit()
    # for sub_name in sorted(well_qc_dict.keys(), key=lambda sub_name: recognized_well_by_file_name(sub_name)):   
    for sub_name in sorted(well_qc_dict.keys(), key=lambda sub_name: sub_name):   
        # print(sub_name)
        # continue
        mis_out = "0"
        if not well_qc_dict[sub_name]["coverage"]:
            mis_out = "-"
        if len(well_qc_dict[sub_name]["mismatch"]) >= 1 and well_qc_dict[sub_name]["coverage"]:
            mis_out = ",".join([str(x) for x in well_qc_dict[sub_name]["mismatch"]])
        # elif len(well_qc_dict[well]['mismatch']) == 1 and well_qc_dict[well]['coverage']:
        #     mis_out = str(well_qc_dict[well]['mismatch'][0])
        all_detective = not 0 in well_qc_dict[sub_name]["sgRNA"]
        well_qc_str = (
            sub_name
            + "\t"
            + str(recognized_well_by_file_name(sub_name))
            + "\t"
            + "\t".join([str(x) for x in well_qc_dict[sub_name]["sgRNA"]])
            + "\t"
            + str(all_detective)
            + "\t"
            + str(well_qc_dict[sub_name]["coverage"])
            + "\t"
            + str(len(well_qc_dict[sub_name]["qc_failed"]))
            + "\t"
            + mis_out
            + "\t"
            + str(len(well_qc_dict[sub_name]["indel_soft"]))
            + "\n"
        )
        out_table_handle.write(well_qc_str)


def pack_ref(output_path: str, output_ref_path: str, output_ref_pack_path: str):
    system(
        f"cp {output_ref_path}/*fa {output_ref_pack_path}/ && tar -cvf {output_path}/ref_pack.tar {output_ref_pack_path} > /dev/null"
    )


def process_pipeline(
    trim_start: int,
    trim_end: int,
    raw_vector_path: str,
    sgRNA_table_path: str,
    input_path: str,
    output_path: str,
) -> None:
    # mkdir and get subplate
    output_fq_path = f"{output_path}/fq/"
    output_ref_path = f"{output_path}/ref/"
    output_ref_pack_path = f"{output_path}/ref_pack/"
    output_bam_path = f"{output_path}/bam/"
    for dir in [output_fq_path, output_ref_path, output_bam_path, output_ref_pack_path]:
        Path(dir).mkdir(exist_ok=1, parents=1)

    # construct reference of each well
    sg_pos_li, vector_li, feature_dict = cut_seq_obtain_pos(raw_vector_path)
    well_ref_dict = generate_ref(sgRNA_table_path, vector_li, output_ref_path)

    pack_ref(output_path, output_ref_path, output_ref_pack_path)
    # sanger file process
    # file_dict = classify_file_by_well(input_path)
    well_qc_dict = defaultdict(lambda: defaultdict(list))
    generate_fq(trim_start, trim_end, well_qc_dict, input_path, output_fq_path)

    # alignment and output
    for fq in Path(output_fq_path).glob("*fq"):
        prefix_li = fq.name.split(".")
        sub_cycle_info = prefix_li[0]
        well = sub_cycle_info.split("_")[1]
        plate = "-".join(sub_cycle_info.split("-")[:2])
        ref_file = f"{output_ref_path}/{plate}_{well}.fa"
        input_fq = str(fq)
        bam_name = fq.name.replace("fq","bam")
        output_bam = f"{output_bam_path}/{bam_name}"
        process_alignment(
            ref_file,
            input_fq,
            output_bam,
            sg_pos_li,
            well_qc_dict,
            feature_dict,
        )
    qc_dict_to_table(
        well_qc_dict, f"{output_path}/res.tsv"
    )




if __name__ == "__main__":
    process_pipeline(
        trim_start=50,
        trim_end=800,
        input_path="/home/wayne/Project/SC/Sanger/0918_test/HA-4-1/",
        sgRNA_table_path="/home/wayne/Project/SC/Sanger/0918_test/HA-4-1/HA-4.xlsx",
        raw_vector_path="/home/wayne/Project/SC/Sanger/pYJA5-4sgRNA.gb",
        output_path="/home/wayne/Project/SC/Sanger/subout/HA-4-1///HA-4_res",
    )
    