import numpy
import argparse
parser=argparse.ArgumentParser(description="Randomly Subsample SNPs from a genepop file")
parser.add_argument("-i", help="Input data file in genepop format", type=str, action='store', required=True)
parser.add_argument("-o", help="Subset output file in genepop format", type=str, action='store', required=True)
parser.add_argument("-n", help="Number of SNPs to subsample", type=int, action='store', required=True)
args=parser.parse_args()

n_snps=args.n
input_data=args.i
output_subset=args.o

snp_file=open(input_data)
subset_snp_file=open(output_subset, 'w')
subset_snp_file.write(snp_file.readline())
loci_names=snp_file.readline().strip().split(",")
random_loci=numpy.random.choice(len(loci_names), n_snps, replace=False)
random_loci_names=[loci_names[i] for i in random_loci]
subset_snp_file.write("\t"+",".join(random_loci_names)+"\n")
for line in snp_file:
    if line.strip()=="POP":
        subset_snp_file.write(line)
    else:
        loci=line.strip().split("\t")
        random_loci_SNPs=[loci[1:][i] for i in random_loci]
        subset_snp_file.write(loci[0]+"\t"+"\t".join(random_loci_SNPs)+"\n")
        print("Individual "+loci[0][0:-1]+" subsampled")
snp_file.close()
subset_snp_file.close()
