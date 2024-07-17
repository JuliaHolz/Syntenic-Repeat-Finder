import argparse 
from Bio import SeqIO
import subprocess
parser = argparse.ArgumentParser(
                    prog='AddBasePairs',
                    description='takes in a sorted bed file of genomic locations, and adds specified number of basepairs to both sides',
                    epilog='Uses pybedtools slop and filter')
parser.add_argument('-t', '--target', action="store", help = "input fasta of repeat locations in target genome")
parser.add_argument('-q', '--query', action="store", help = "input fasta of repeat locations in query genome (mapped over from target)")
parser.add_argument('-o', '--output', action="store", help = "output file, where to output the alignment")

args = parser.parse_args()


target_fasta = SeqIO.parse(args.target, "fasta")
query_fasta =  SeqIO.parse(args.query, "fasta")

for target_record in target_fasta:
    query_record = next(query_fasta)
    if query_record.id != target_record.id:
        raise Exception("Error: target and query names do not match")
    else:
        subprocess.run("")