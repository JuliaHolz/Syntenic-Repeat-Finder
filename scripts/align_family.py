import subprocess
import argparse 
import time

parser = argparse.ArgumentParser(
                    prog='AddBasePairs',
                    description='takes in a sorted bed file of genomic locations, and filters out those without unique sequence on one or more sides, outputting a new, filtered bed file.',
                    epilog='Uses pybedtools  filter')
parser.add_argument('-i', '--input', action="store", help = "bed file containing the repeats in the target genome")
parser.add_argument('-m', '--mapped', action="store", help = "file containing the repeats mapped to the query genome")
parser.add_argument('-t', '--targetgenome', action="store", help = "fasta file containing the target genome (which -i bed file refers to)")
parser.add_argument('-q', '--querygenome', action="store", help = "fasta file containing the query genome (which -m bed file refers to)")
parser.add_argument('-o', '--output', action="store", help = "location where outputs for this genome/set of repeats should be output (alignments, original_fasta and mapped_fasta directories should already exist in this location)")
parser.add_argument('-f', '--family', action="store", help = "family name")

args = parser.parse_args()
target_genome_name = (args.targetgenome.split("/")[-1]).split(".")[0]
query_genome_name = (args.querygenome.split("/")[-1]).split(".")[0]
start_fasta = time.time()
#get fasta files for the bed
print("converting target repeats to fasta")
target_output_file = args.output + "/original_fasta/" + args.family + ".fasta"
subprocess.run(["bedtools", "getfasta", "-fi", args.targetgenome, "-bed", args.input, "-fo", target_output_file, "-name"])
query_output_file = args.output + "/mapped_fasta/" + args.family + ".fasta"
print("converting query repeats to fasta")
subprocess.run(["bedtools", "getfasta", "-fi", args.querygenome, "-bed", args.mapped, "-fo", query_output_file, "-name"])
end_fasta = time.time()
alignment_output_file = args.output + "/alignments/" + args.family + "/alignments.caf"

with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "w") as outfile:
    outfile.write("created fasta from bed in: " + str(end_fasta-start_fasta) + " seconds \n")


start_align = time.time()
print("aligning repeats")
target_file =open(target_output_file,"r")
query_file =open(query_file,"r")
alignment_file = open(alignment_output_file, "w")

opened = False
for target_line in target_file:
    query_line = query_file.readline()
    if(target_line[0] == ">") :
        if(opened):
            target_temp.close()
            query_temp.close()
            subprocess.run( ["/usr/local/RepeatModeler/util/align.pl", "-force", "-crossmatch", "-caf", target_fasta, query_fasta], stdout=alignment_file)
            subprocess.run( ["rm", target_temp.name])
            subprocess.run( ["rm", query_temp.name])
        opened = True
        filename = (target_line[1:].rstrip())
        query_filename = (query_line[1:].rstrip())
        if(query_filename != filename):
            raise Exception("mismatch between target and query entry names {} and {}".format(filename,query_filename))
        filename = filename.replace("/", "%")
        filename = filename.replace(";", "#")
        filename = filename.replace("(", "@")
        filename = filename.replace(")", "@")
        target_line = ">" + target_genome_name + ":" + filename.split(":")[0] + "\n"
        query_line = ">" + target_genome_name + ":" + filename.split(":")[0] + "\n"
        target_temp=open(args.output + "/" + target_genome_name + "%s.fa" % filename, "w")
        query_temp=open(args.output + "/" + query_genome_name + "%s.fa" % filename, "w")
    target_temp.write(target_line)
target_temp.close()
query_temp.close()
subprocess.run( ["rm", target_temp.name])
subprocess.run( ["rm", query_temp.name])
alignment_file.close()

end_align = time.time()

with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "a") as outfile:
    outfile.write("aligned " +  str(len(repeat_names)) +  " repeats in " + str(end_align-start_align) + " seconds \n")
