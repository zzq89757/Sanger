from Bio import SeqIO

# seq = SeqIO.parse("/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1","abi")

# print(dir(seq))

for seq in SeqIO.parse("/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1","abi"):
        seq_id, sequence, qual_array = seq.name, seq.seq, seq.letter_annotations['phred_quality']
        abif_raw = seq.annotations["abif_raw"]
        data_g = list(abif_raw["DATA9"])
        data_a = list(abif_raw["DATA10"])
        data_t = list(abif_raw["DATA11"])
        data_c = list(abif_raw["DATA12"])

        print(data_g)