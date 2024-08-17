import subprocess
import argparse 
import time

parser = argparse.ArgumentParser(
                    prog='download_consensus',
                    description='downloads consensus from Dfam',
                    epilog='using Dfam REST API')
parser.add_argument('-f', '--fasta', action="store", help = "family name")
parser.add_argument('-o', '--output', action="store", help = "output file")
args = parser.parse_args()

with open(args.fasta, "r") as input_fasta:
    with open(args.output, "w") as output_fasta:
        for line in input_fasta:
            if line[0] == ">":
                split = line.split(";")
                second_section = split[3].split(":")[0]
                id = ">" + split[2] + ";" +  second_section + ";" + split[3].split(":")[2]
                output_fasta.write(id + "\n")
            else:
                output_fasta.write(line)
