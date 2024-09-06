from math import ceil
from Bio  import SeqIO, SeqRecord
def fq_from_abi(file_path: str) -> None:
    '''读取ab1文件并将序列信息存入fq格式字符串'''
    fq_str = ""
    for seq in SeqIO.parse(file_path, "abi"):
        seq_id = seq.name
        # base signal
        abif_raw = seq.annotations["abif_raw"]
        peak_loc_li = abif_raw['PLOC2']
        # print(len(peak_loc_li))
        print(len(seq))
        a_trace = list(abif_raw["DATA9"])
        print(len(a_trace))
        g_trace = list(abif_raw["DATA10"])
        t_trace = list(abif_raw["DATA11"])
        c_trace = list(abif_raw["DATA12"])
        trace_length = len(a_trace)
        loc = 1
        loc_idx = 0
        
        # vectorbee method
        start_pos: int = 0
        next_base_pos: int = peak_loc_li[1]
        end_pos: int = 0
        
        def set_end_pos(start_pos: int, next_base_pos: int, end_pos: int, trace_len: int) -> int:
            if next_base_pos:
                end_pos = start_pos + ceil((next_base_pos - start_pos) / 2)
            else:
                end_pos = trace_len
            return end_pos
        end_pos = set_end_pos(start_pos, next_base_pos, end_pos)
        print(len(peak_loc_li))
        for i in range(len(peak_loc_li)):
            print(f"loc: {i}")
            for j in range(start_pos,end_pos + 3):
                # if j == len(a_trace):break
                a, g, c, t = a_trace[j], g_trace[j], c_trace[j], t_trace[j]
                print(f"{j}-A:{a}\tT:{t}\tC:{c}\tG:{g}")
            loc += 1
            # next_base_pos = None
            if i != len(peak_loc_li) - 1:
                start_pos = end_pos + 3
                next_base_pos = peak_loc_li[i + 2] if i != len(peak_loc_li) - 2 else None
                end_pos = set_end_pos(start_pos, next_base_pos, end_pos)
                print(end_pos)
            
            
        # for i in range(len(peak_loc_li) -1 ):
        #     print(f"loc: {loc}")
        #     start_idx = peak_loc_li[i] - 1 
        #     end_idx = peak_loc_li[i + 1] - 1
        #     for j in range(start_idx,end_idx):
        #         a, g, c, t = a_trace[j], g_trace[j], c_trace[j], t_trace[j]
        #         print(f"{a}\t{t}\t{c}\t{g}")
        #     loc += 1
        # for a,g,c,t in zip(a_trace,g_trace,c_trace,t_trace):
        #     if loc == peak_loc_li[loc_idx]:
        #         print("\n")
        #         loc_idx += 1
        #     print(f"{a}\t{t}\t{c}\t{g}")
            # loc += 1
        # signal qc
        # for 
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