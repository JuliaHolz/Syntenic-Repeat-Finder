import subprocess
import argparse 
import time

parser = argparse.ArgumentParser(
                    prog='AlignFamily',
                    description='aligns all instances of a family using crossmatch, and outputs them in caf/yaCaf format',
                    epilog='outputs target and query fastas with sequences of repeats, as alignments.caf (which contains all the alignments) ')
parser.add_argument('-i', '--input', action="store", help = "bed file containing the repeats in the target genome")
parser.add_argument('-m', '--mapped', action="store", help = "file containing the repeats mapped to the query genome")
parser.add_argument('-t', '--targetgenome', action="store", help = "fasta file containing the target genome (which -i bed file refers to)")
parser.add_argument('-q', '--querygenome', action="store", help = "fasta file containing the query genome (which -m bed file refers to)")
parser.add_argument('-o', '--output', action="store", help = "location where outputs for this genome/set of repeats should be output (alignments, target_fasta and query_fasta directories should already exist in this location) and alignments_summary.txt")
parser.add_argument('-f', '--family', action="store", help = "family name")

args = parser.parse_args()
target_genome_name = (args.targetgenome.split("/")[-1]).split(".")[0]
query_genome_name = (args.querygenome.split("/")[-1]).split(".")[0]
start_fasta = time.time()
#get fasta files for the bed
print("converting target repeats to fasta")
target_output_file = args.output + "/target_fasta/" + args.family + ".fasta"
subprocess.run(["bedtools", "getfasta", "-fi", args.targetgenome, "-bed", args.input, "-fo", target_output_file, "-name"])
query_output_file = args.output + "/query_fasta/" + args.family + ".fasta"
print("converting query repeats to fasta")
subprocess.run(["bedtools", "getfasta", "-fi", args.querygenome, "-bed", args.mapped, "-fo", query_output_file, "-name"])
end_fasta = time.time()
alignment_output_file = args.output + "/alignments/" + args.family + "/alignments.caf"
error_file_name = args.output + "/alignments/" + args.family + "/errors.txt"

with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "w") as outfile:
    outfile.write("created fasta from bed in: " + str(end_fasta-start_fasta) + " seconds \n")

num_repeats = 0
start_align = time.time()
print("aligning repeats")
target_file =open(target_output_file,"r")
query_file =open(query_output_file,"r")
alignment_file = open(alignment_output_file, "w")
error_file = open(error_file_name, "w")
opened = False
def align_repeats(target_filename,query_filename, alignment_file):
    alignment_result = subprocess.run( ["/usr/local/RepeatModeler/util/align.pl", "-force", "-crossmatch", "-caf", target_filename, query_filename], stdout=subprocess.PIPE, stderr=error_file)
    alignment_lines = alignment_result.stdout.decode("utf-8").split("\n")
    num_alignments_expected = int([line for line in alignment_lines if line[0:13] == "#  alignments"][0].split(":")[1])
    alignments = [line for line in alignment_lines if line != '' and line[0]!="#"]
    assert(len(alignments) == num_alignments_expected)
    alignment_file.write(">" + str(num_alignments_expected) + " alignments for:" + filename.split(":")[0] + "\n")
    for alignment in alignments:
        alignment_file.write(alignment + "\n")

for target_line in target_file:
    query_line = query_file.readline()
    if(target_line[0] == ">") :
        num_repeats +=1
        if(opened):
            log_file_name = args.output + "/" + target_genome_name + "%s.fa.log" % filename
            target_temp.close()
            query_temp.close()
            align_repeats(target_filename,query_filename,alignment_file)

            subprocess.run( ["rm", "-f", target_temp.name])
            subprocess.run( ["rm", "-f", query_temp.name])
            subprocess.run( ["rm", "-f", log_file_name])
        opened = True
        filename = (target_line[1:].rstrip())
        query_filename = (query_line[1:].rstrip())
        if(query_filename.split(":")[0] != filename.split(":")[0]):
            raise Exception("mismatch between target and query entry names {} and {}".format(filename,query_filename))
        filename = filename.replace("/", "%")
        filename = filename.replace(";", "#")
        filename = filename.replace("(", "@")
        filename = filename.replace(")", "@")
        query_filename = query_filename.replace("/", "%")
        query_filename = query_filename.replace(";", "#")
        query_filename = query_filename.replace("(", "@")
        query_filename = query_filename.replace(")", "@")
        target_line = ">" + target_genome_name + ":" + filename.split(":")[0] + "\n"
        query_line = ">" + query_genome_name + ":" + filename.split(":")[0] + "\n"
        target_filename = args.output + "/" + target_genome_name + "%s.fa" % filename
        query_filename = args.output + "/" + query_genome_name + "%s.fa" % query_filename
        target_temp=open(target_filename, "w")
        query_temp=open(query_filename, "w")
    target_temp.write(target_line)
    query_temp.write(query_line)
    
log_file_name = args.output + "/" + target_genome_name + "%s.fa.log" % filename
target_temp.close()
query_temp.close()

align_repeats(target_filename,query_filename,alignment_file)

subprocess.run( ["rm", "-f", target_temp.name])
subprocess.run( ["rm", "-f", query_temp.name])
subprocess.run( ["rm", "-f", log_file_name])

error_file.close()

alignment_file.close()

end_align = time.time()

with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "a") as outfile:
    outfile.write("aligned " +  str(num_repeats) +  " repeats in " + str(end_align-start_align) + " seconds \n")
