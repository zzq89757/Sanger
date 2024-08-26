from Bio import SeqIO

seq = SeqIO.parse("/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1","abi")

# print(dir(seq))

for s in seq:
  print(s.name)
  print(s.seq)
  # phred_quality
  print(s.letter_annotations['phred_quality'])