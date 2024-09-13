from Bio import SeqIO

def parse_gb_file_to_dict(gb_file: str) -> dict:
    """
    解析 GenBank 文件,将 feature 的起始和结束位置及其 /label 存入字典。

    :param gb_file: GenBank 文件路径
    :return: 包含位置和 /label 的字典
    """
    feature_dict = {}
    for record in SeqIO.parse(gb_file, "genbank"):
        # 遍历所有的 features
        for feature in record.features:
            start = int(feature.location.start)
            end = int(feature.location.end)
            
            # 如果 feature 有 /label,则存入字典
            if "label" in feature.qualifiers:
                label = feature.qualifiers['label'][0]
                # 将位置范围作为键,label 作为值存入字典
                feature_dict[(start, end)] = label
    return feature_dict

def find_label_by_position(feature_dict: dict, position: int) -> str:
    """
    根据输入的位置从字典中查找对应的 /label。

    :param feature_dict: feature 位置信息和 /label 的字典
    :param position: 输入的位置（从 0 开始）
    :return: 对应位置的 label,若未找到则返回 "No label found."
    """
    for (start, end), label in feature_dict.items():
        if start <= position <= end:
            return label
    return "No label found."



# 测试函数
gb_file = "/home/wayne/Project/SC/Sanger/pYJA5-4sgRNA.gb"  # 替换为你的 GenBank 文件路径
feature_dict = parse_gb_file_to_dict(gb_file)
position = 5200  # 输入的查找位置
label = find_label_by_position(feature_dict, position)
print(f"Position {position} corresponds to label: {label}")
