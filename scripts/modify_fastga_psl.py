import argparse 

parser = argparse.ArgumentParser(
                    prog='RemoveDescriptionsFromPSL',
                    description='removes verbose descriptions from a psl file output by FASTGA')
parser.add_argument('-i', '--input', action="store", help = "psl file")
parser.add_argument('-o', '--output', action="store", help = "location to output the new psl file to")

args = parser.parse_args()
inpath = args.input
outpath = args.output
count =1
with open(inpath, "r") as pslfile:
    with open(outpath, "w") as outfile:
        for line in pslfile:
            count+=1
            line=line.split("\t")
            #print(list(zip(line, range(0, len(line)))))
            line[9] = line[9].split()[0]
            line[13] = line[13].split()[0]
            #line[10] = str(int(line[11])-1)
            #line[12] = str(int(line[12])-1)
            #line[14] = str(int(line[15])-1)
            #line[16] = str(int(line[16])-1)
            #subtract 1 from the query and target starts because FASTGA does not output them zero-indexed
            qstarts = [str(int(s)-1) for s in line[19].split(",") if (len(s)>0 and not s.isspace())]
            line[19] = ",".join(qstarts)
            tstarts = [str(int(s)-1) for s in line[20].split(",") if (len(s)>0 and not s.isspace())]
            line[20] = ",".join(tstarts)
            outfile.write("\t".join(line))
            #if count>=5: 
            #    break
 