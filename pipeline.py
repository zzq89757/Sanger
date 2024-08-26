from Bio import SeqIO


def parse_abi_file(file_path):
    for seq in SeqIO.parse(file_path,"abi"):
        print(seq.name)
        print(seq.seq)


def main():
    input_file = "/home/wayne/Project/SC/Sanger/B103-(T3389)pUp-pDown-flank-R.ab1"
    parse_abi_file(input_file)
    

if __name__ == "__main__":
    main()