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
start_fasta = time.time()
#get fasta files for the bed
print("converting target repeats to fasta")
target_output_file = args.output + "/original_fasta/" + args.family + ".fasta"
subprocess.run(["bedtools", "getfasta", "-fi", args.targetgenome, "-bed", args.input, "-fo", target_output_file, "-name"])
query_output_file = args.output + "/mapped_fasta/" + args.family + ".fasta"
print("converting query repeats to fasta")
subprocess.run(["bedtools", "getfasta", "-fi", args.querygenome, "-bed", args.mapped, "-fo", query_output_file, "-name"])
end_fasta = time.time()

with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "a") as outfile:
    outfile.write("created fasta from bed in: " + str(end_fasta-start_fasta) + " seconds \n")

subprocess.run(["mkdir", "-p", args.output + "/alignments/" + args.family])
target_repeat_folder =  args.output + "/alignments/" + args.family + "/target_repeats/"
subprocess.run(["mkdir", "-p", target_repeat_folder])
query_repeat_folder = args.output + "/alignments/" + args.family + "/query_repeats/"
subprocess.run(["mkdir", "-p",query_repeat_folder])
aligned_output_folder = args.output + "/alignments/" + args.family + "/aligned/"
subprocess.run(["mkdir", "-p", aligned_output_folder])

repeat_names = []
start_split = time.time()
print("splitting target repeats fasta file into individual fastas")
f=open(target_output_file,"r")
opened = False
for line in f :
    if(line[0] == ">") :
        if(opened):
            of.close()
        opened = True
        filename = (line[1:].split(":")[0].rstrip())
        filename = filename.replace("/", "%")
        filename = filename.replace(";", "#")
        filename = filename.replace("(", "@")
        filename = filename.replace(")", "@")
        line = ">" + filename + "\n"
        repeat_names.append(filename)
        of=open(target_repeat_folder + "%s.fa" % filename, "w")
    of.write(line)
of.close()

print("splitting query repeats fasta file into individual fastas")
f=open(query_output_file,"r")
opened = False
for line in f :
    if(line[0] == ">") :
        if(opened):
            of.close()
        opened = True
        filename = (line[1:].split(":")[0].rstrip())
        filename = filename.replace("/", "%")
        filename = filename.replace(";", "#")
        filename = filename.replace("(", "@")
        filename = filename.replace(")", "@")
        line = ">" + filename + "\n"
        of=open(query_repeat_folder + "%s.fa" % filename, "w")
    of.write(line)
of.close()
end_split = time.time()
with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "a") as outfile:
    outfile.write("split fastas into individual repeats in: " + str(end_split-start_split) + " seconds \n")
start_align = time.time()
print("aligning repeats")
count = 0
for repeatName in repeat_names:
    target_fasta = target_repeat_folder + "%s.fa" % repeatName
    query_fasta = query_repeat_folder + "%s.fa" % repeatName
    output_location = aligned_output_folder + repeatName + ".txt"
    with open(output_location, "w") as outfile:
        proc = subprocess.run( ["/usr/local/RepeatModeler/util/align.pl", "-force", "-crossmatch", "-a", target_fasta, query_fasta], stdout=outfile)
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.nhr"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.nin"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.njs"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.nog"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.nos"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.not"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.nsq"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.ntf"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.nto"])
        #subprocess.run( ["rm",  args.output + "/alignments/" + args.family + "/query_repeats/"+repeatName+".fa.ndb"])
#create the summary file
end_align = time.time()

with open(args.output + "/alignments/" + args.family + "/alignment_summary.txt", "a") as outfile:
    outfile.write("aligned " +  str(len(repeat_names)) +  " repeats in " + str(end_align-start_align) + " seconds \n")
