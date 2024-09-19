from pysam import AlignmentFile, AlignedSegment
from collections import defaultdict

handle = AlignmentFile("/home/wayne/Project/SC/Sanger/HA-1_res/bam/HA-1-4-B07_86.bam")
def found_mismatch_by_md_tag(
    well: int,
    aln: AlignedSegment,
    well_qc_dict: defaultdict,
    detective_end: int,
    feature_dict: dict,
) -> bool:
    md_str = aln.get_tag("MD")
    print(md_str)
    forward_length = ''
    abs_pos = aln.reference_start
    forward_sum = 0
    for i in md_str:
        # obtain ascii code 
        if ord(i) <= 65: # 0~9
            forward_length += i
        else:
            # current snp position and label
            abs_pos += int(forward_length) + 1
            forward_sum += int(forward_length)
            forward_length = ''
            # component = find_label_by_position(feature_dict, abs_pos)
            component = 'MVC'
            # 
            print(forward_sum)
            mismatch_str = f"{abs_pos}<{component}>:{i}->{aln.query_alignment_sequence[forward_sum]}"
            if abs_pos < detective_end and mismatch_str not in well_qc_dict[well]["mismatch"]:
                # stop codon check
                # print(aln.query_name)
                # stop_codon_check(
                #     aln.query_alignment_sequence[forward_sum - 2 : forward_sum + 3],
                #     well_qc_dict,
                # )
                well_qc_dict[well]["mismatch"].append(
                    f"{abs_pos}<{component}>:{i}->{aln.query_alignment_sequence[forward_sum]}"
                )
            forward_sum += 1
well_qc_dict = defaultdict(lambda: defaultdict(list))
feature_dict = dict()
for aln in handle:
    cs = aln.get_cigar_stats()
    cs = aln.get_tag("MD")
    found_mismatch_by_md_tag(1,aln,well_qc_dict,9000,feature_dict)

print(well_qc_dict)
    
    
    
