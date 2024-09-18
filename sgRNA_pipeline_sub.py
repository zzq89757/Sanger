from collections import defaultdict
from os import system
from pathlib import Path
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment
import pandas as pd


def cut_seq_obtain_pos(seq_path: str) -> list:
    """根据原始载体构建sgRNA位置列表以及载体不含sgRNA的切片列表和元件位置字典"""
    # for v in SeqIO.parse(vector_seq, 'fasta'):
    #     vector_seq = str(v.seq)
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
    subplate, raw_col = file_name.split("-")[2:4]
    raw_num = raw_col[0]
    col_num = int(raw_col[1:])
    raw_num = ord(raw_num) - 65
    sub_well_idx = raw_num * 12 + int(col_num) - 1
    # print(subplate)
    print(file_name.split("-")[2:4])
    well = sub_dict[subplate][sub_well_idx]
    return well
for a in ['HA-1-4-G12-P3_G12', 'HA-1-2-E07-P3_E07', 'HA-1-4-A09-P1_A09', 'HA-1-3-D01-P3_D01', 'HA-1-3-A12-P1_A12', 'HA-1-2-G07-P3_G07', 'HA-1-4-F04-P3_F04', 'HA-1-4-F10-P1_F01', 'HA-1-4-G10-P3_G10', 'HA-1-2-E01-P3_E01', 'HA-1-4-G10-P1_G10', 'HA-1-2-C01-P2_C01', 'HA-1-3-E12-P3_E12', 'HA-1-3-E12-P2_E12', 'HA-1-3-F01-P2_F01', 'HA-1-4-F02-P3_F02', 'HA-1-3-F08-P2_F08', 'HA-1-3-C12-P3_C12', 'HA-1-3-C07-P1_C07', 'HA-1-4-C07-P1_C07', 'HA-1-4-D07-P1_D07', 'HA-1-3-H06-P3_H06', 'HA-1-2-E12-P3_E12', 'HA-1-4-H07-P1_H07', 'HA-1-4-B02-P1_B02', 'HA-1-2-F05-P3_F05', 'HA-1-3-H10_379', 'HA-1-2-D07_158', 'HA-1-4-G12_336', 'HA-1-2-D08_160', 'HA-1-4-F01_266', 'HA-1-3-D11_189', 'HA-1-4-G07_326', 'HA-1-3-B02_75', 'HA-1-4-F10_284', 'HA-1-3-F08_279', 'HA-1-2-H09_354', 'HA-1-2-H07_350', 'HA-1-4-A07_38', 'HA-1-3-H07_373', 'HA-1-2-E05_202', 'HA-1-4-H04_368', 'HA-1-2-B11_70', 'HA-1-4-D07_182', 'HA-1-3-C10_139', 'HA-1-4-B11_94', 'HA-1-2-H05_346', 'HA-1-2-H11_358', 'HA-1-3-F07_277', 'HA-1-3-F04_271', 'HA-1-2-G04_296', 'HA-1-3-H06_371', 'HA-1-2-F02_244', 'HA-1-4-B04_80', 'HA-1-4-A01_26', 'HA-1-4-H11_382', 'HA-1-4-F08_280', 'HA-1-3-A01_25', 'HA-1-2-H06_348', 'HA-1-4-D10_188', 'HA-1-2-G10_308', 'HA-1-2-D12_168', 'HA-1-4-H12_384', 'HA-1-3-D01_169', 'HA-1-2-G03_294', 'HA-1-4-D02_172', 'HA-1-2-B01_50', 'HA-1-4-B12_96', 'HA-1-3-D06_179', 'HA-1-2-B06_60', 'HA-1-4-D09_186', 'HA-1-2-B07_62', 'HA-1-3-B09_89', 'HA-1-3-E10_235', 'HA-1-4-H10_380', 'HA-1-2-G08_304', 'HA-1-4-F03_270', 'HA-1-4-E05_226', 'HA-1-4-F12_288', 'HA-1-2-F01_242', 'HA-1-2-D09_162', 'HA-1-3-G03_317', 'HA-1-4-C07_134', 'HA-1-4-D12_192', 'HA-1-2-B04_56', 'HA-1-4-F09_282', 'HA-1-4-D03_174', 'HA-1-2-H08_352', 'HA-1-4-E02_220', 'HA-1-4-G04_320', 'HA-1-3-C02_123', 'HA-1-2-E01_194', 'HA-1-2-C09_114', 'HA-1-4-H05_370', 'HA-1-4-D06_180', 'HA-1-3-A09_41', 'HA-1-4-C04_128', 'HA-1-2-E07_206', 'HA-1-2-C12_120', 'HA-1-2-E06_204', 'HA-1-3-H12_383', 'HA-1-2-F06_252', 'HA-1-3-F02_267', 'HA-1-4-B06_84', 'HA-1-3-B08_87', 'HA-1-3-F11_285', 'HA-1-2-H01_338', 'HA-1-4-D05_178', 'HA-1-2-G12_312', 'HA-1-2-G07_302', 'HA-1-4-G08_328', 'HA-1-4-H01_362', 'HA-1-3-D02_171', 'HA-1-3-G01_313', 'HA-1-3-C07_133', 'HA-1-4-G02_316', 'HA-1-2-C08_112', 'HA-1-4-G01_314', 'HA-1-4-F04_272', 'HA-1-2-C05_106', 'HA-1-4-A02_28', 'HA-1-2-C01_98', 'HA-1-3-C04_127', 'HA-1-2-F05_250', 'HA-1-4-G11_334', 'HA-1-2-E04_200', 'HA-1-3-E01_217', 'HA-1-4-A03_30', 'HA-1-2-B03_54', 'HA-1-3-F01_265', 'HA-1-2-C02_100', 'HA-1-3-B04_79', 'HA-1-2-G05_298', 'HA-1-3-D10_187', 'HA-1-4-D11_190', 'HA-1-2-C04_104', 'HA-1-4-E11_238', 'HA-1-2-F07_254', 'HA-1-3-B12_95', 'HA-1-4-E12_240', 'HA-1-4-C06_132', 'HA-1-3-H09_377', 'HA-1-4-C12_144', 'HA-1-3-D08_183', 'HA-1-3-C11_141', 'HA-1-4-E10_236', 'HA-1-3-A02_27', 'HA-1-2-D06_156', 'HA-1-3-G07_325', 'HA-1-4-H07_374', 'HA-1-2-D05_154', 'HA-1-2-B05_58', 'HA-1-3-E06_227', 'HA-1-2-F11_262', 'HA-1-4-E07_230', 'HA-1-3-H05_369', 'HA-1-3-E08_231', 'HA-1-4-G09_330', 'HA-1-3-G10_331', 'HA-1-2-G09_306', 'HA-1-3-A07_37', 'HA-1-4-H03_366', 'HA-1-2-G11_310', 'HA-1-3-G09_329', 'HA-1-2-H02_340', 'HA-1-2-C03_102', 'HA-1-3-F09_281', 'HA-1-3-B10_91', 'HA-1-4-A06_36', 'HA-1-3-A10_43', 'HA-1-3-E11_237', 'HA-1-2-F10_260', 'HA-1-2-D11_166', 'HA-1-3-G12_335', 'HA-1-2-H04_344', 'HA-1-3-H08_375', 'HA-1-4-G05_322', 'HA-1-2-D03_150', 'HA-1-2-D02_148', 'HA-1-4-A04_32', 'HA-1-4-B09_90', 'HA-1-3-B05_81', 'HA-1-2-F09_258', 'HA-1-3-B11_93', 'HA-1-3-F06_275', 'HA-1-2-E10_212', 'HA-1-4-H09_378', 'HA-1-3-D12_191', 'HA-1-3-E12_239', 'HA-1-4-D08_184', 'HA-1-4-G06_324', 'HA-1-2-H10_356', 'HA-1-2-E02_196', 'HA-1-4-D04_176', 'HA-1-2-F04_248', 'HA-1-4-E04_224', 'HA-1-3-C09_137', 'HA-1-4-A05_34', 'HA-1-3-C05_129', 'HA-1-2-B02_52', 'HA-1-4-A12_48', 'HA-1-3-D04_175', 'HA-1-3-E09_233', 'HA-1-2-C07_110', 'HA-1-3-H04_367', 'HA-1-2-E11_214', 'HA-1-3-G06_323', 'HA-1-3-F10_283', 'HA-1-4-B02_76', 'HA-1-3-C01_121', 'HA-1-2-B08_64', 'HA-1-3-H02_363', 'HA-1-2-B09_66', 'HA-1-4-B10_92', 'HA-1-4-F02_268', 'HA-1-4-A10_44', 'HA-1-4-F06_276', 'HA-1-3-A06_35', 'HA-1-2-B10_68', 'HA-1-2-F08_256', 'HA-1-3-D05_177', 'HA-1-3-B06_83', 'HA-1-2-A09_18', 'HA-1-4-B01_74', 'HA-1-4-H08_376', 'HA-1-3-E05_225', 'HA-1-2-G01_290', 'HA-1-2-H03_342', 'HA-1-4-B03_78', 'HA-1-2-E09_210', 'HA-1-4-C01_122', 'HA-1-3-A12_47', 'HA-1-2-C06_108', 'HA-1-4-B05_82', 'HA-1-3-E07_229', 'HA-1-3-G11_333', 'HA-1-2-E03_198', 'HA-1-3-E02_219', 'HA-1-3-G02_315', 'HA-1-2-F03_246', 'HA-1-4-E06_228', 'HA-1-3-A08_39', 'HA-1-4-H06_372', 'HA-1-4-F07_278', 'HA-1-2-E12_216', 'HA-1-3-E04_223', 'HA-1-4-G10_332', 'HA-1-3-A11_45', 'HA-1-3-B03_77', 'HA-1-2-C10_116', 'HA-1-3-G08_327', 'HA-1-3-H03_365', 'HA-1-4-A09_42', 'HA-1-3-A03_29', 'HA-1-2-D10_164', 'HA-1-3-C12_143', 'HA-1-4-C10_140', 'HA-1-3-B01_73', 'HA-1-4-C03_126', 'HA-1-4-E08_232', 'HA-1-4-C11_142', 'HA-1-2-D01_146', 'HA-1-2-G06_300', 'HA-1-4-F05_274', 'HA-1-3-A05_33', 'HA-1-4-E03_222', 'HA-1-4-E09_234', 'HA-1-4-C09_138', 'HA-1-2-E08_208', 'HA-1-4-B08_88', 'HA-1-3-G04_319', 'HA-1-4-A08_40', 'HA-1-4-C08_136', 'HA-1-3-H11_381', 'HA-1-2-G02_292', 'HA-1-2-C11_118', 'HA-1-2-B12_72', 'HA-1-4-D01_170', 'HA-1-3-A04_31', 'HA-1-4-C02_124', 'HA-1-3-H01_361', 'HA-1-4-G03_318', 'HA-1-3-E03_221', 'HA-1-3-F05_273', 'HA-1-2-D04_152', 'HA-1-4-B07_86']:
    # print(a)
    ...
well = recognized_well_by_file_name('HA-1-3-H10_379')
print(well)
exit()


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
            well_qc_dict[seq_id]["qc_failed"] = [1]
            print(f"{seq_id} no signal !!!")
            return ""
        # trimmed qual control
        if not trimmed_qual_qc(trimmed_qual):
            well_qc_dict[seq_id]["qc_failed"] = [1]
            print(f"{seq_id} qc unpass !!!")
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
        sub_name = file.name.split("-P")[0]
        well = recognized_well_by_file_name(sub_name)
        output_handle = open(f"{output_fq_path}/{sub_name}_{well}.fq", "a")
        fq_str = fq_from_abi(trim_start, trim_end, well_qc_dict, file)
        output_handle.write(fq_str)


def stop_codon_check(sub_name: str, segment: str, well_qc_dict: defaultdict) -> None:
    """检测突变后的碱基是否导致终止密码子的产生"""
    stop_codon_li = ["UAA", "UAG", "UGA"]
    # stop codon detective by ORF start and SNP position

    # if not in cds, continue

    # in cds, get ORF and codon after mutation
    for i in range(3):
        codon = segment[i : i + 3]
        # in stop codon li and in cds
        if codon in stop_codon_li:
            well_qc_dict[sub_name]["stop_codon_from_snp"] = [1]


def mismatch_check(
    sub_name: str,
    aln: AlignedSegment,
    well_qc_dict: defaultdict,
    detective_end: int,
    feature_dict: dict,
) -> bool:
    """检查错配数目,等于1报点,大于1报错"""
    all_snp_count = sum(aln.get_cigar_stats()[0][1:])
    if all_snp_count > 1:
        print(f"{aln.query_name} has mismatch or softclip !!!")
        well_qc_dict[sub_name]["indel_soft"] += [1]
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
        if pos < detective_end and mismatch_str not in well_qc_dict[sub_name]["mismatch"]:
            # stop codon check
            # print(aln.query_name)
            stop_codon_check(
                sub_name,
                aln.query_alignment_sequence[forward_len - 2 : forward_len + 3],
                well_qc_dict,
            )
            well_qc_dict[sub_name]["mismatch"].append(
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
    sub_name = Path(output_bam).name.split(".")[0]
    well_qc_dict[sub_name]["sgRNA"] = [0] * 4
    well_qc_dict[sub_name]["coverage"] = 0
    aln_region_li = []
    for aln in AlignmentFile(output_bam, "r", threads=16):
        # skip supplementary and secondary
        if aln.is_supplementary or aln.is_secondary:
            continue
        # mismatch check
        if mismatch_check(sub_name, aln, well_qc_dict, sg_pos_li[-1] + 100, feature_dict):
            continue
        # sgRNA detective
        cover_idx_li: list[int] = sgRNA_detective(aln, sg_pos_li)
        for i in cover_idx_li:
            well_qc_dict[sub_name]["sgRNA"][i] = 1
        aln_region_li.append([aln.reference_start, aln.reference_end])
    # coverage check after all fastq mapped to ref
    if coverage_check(5914, sg_pos_li[-1] + 150, aln_region_li):
        well_qc_dict[sub_name]["coverage"] = 1


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
    print(well_qc_dict.keys())
    for i in well_qc_dict.keys():
        print(i)
        recognized_well_by_file_name(i)
    exit()
    for sub_name in sorted(well_qc_dict.keys(), key=lambda sub_name: recognized_well_by_file_name(sub_name)):   
        print(sub_name)
        continue
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
        sub_name, well = fq.name.split(".")[0].split("_")
        plate = "-".join(sub_name.split("-")[:2]).upper()
        ref_file = f"{output_ref_path}/{plate}_{well}.fa"
        input_fq = str(fq)
        output_bam = f"{output_bam_path}/{sub_name}_{well}.bam"
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
        input_path="/home/wayne/Project/SC/Sanger/0911/HA-1_raw/",
        sgRNA_table_path="/home/wayne/Project/SC/Sanger/HA-1.xlsx",
        raw_vector_path="/home/wayne/Project/SC/Sanger/pYJA5-4sgRNA.gb",
        output_path="/home/wayne/Project/SC/Sanger/subout//HA-1_res",
    )
    