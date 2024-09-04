from collections import defaultdict, deque
from pathlib import Path


# 提取同一个well中的多个ab1文件中的序列信息 保存为fq后 bwa比对 处理bam文件
def classify_file_by_well(file_path:Path) -> defaultdict:
    '''输入下机数据路径,将其分配至对应的well'''
    ...



def extract_data() -> defaultdict:
    '''提取每个well中的序列信息'''
    ...


def generate_fq():
    '''根据提取的信息生成fq文件(质量值默认为F)'''
    ...
    
    
def detective_alignment_result():
    '''检测是否存在异常的结果'''
    ...
    
    
def process_alignment():
    '''执行bwa命令后处理bam文件'''
    ...
    