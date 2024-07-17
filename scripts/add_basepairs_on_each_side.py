import argparse 
from pybedtools import BedTool
#This script takes in a  bed file of repeat locations

#options when running this script
#-b, how many basepairs to add on each side of each feature 
#-i, input file -- this should be the .bed file of repeat locations
#-g, genome file -- specifies where the genome file (as described in https://bedtools.readthedocs.io/en/latest/content/tools/slop.html) is located
#-o, output file -- where to output the updated .bed file
parser = argparse.ArgumentParser(
                    prog='AddBasePairs',
                    description='takes in a sorted bed file of genomic locations, and adds specified number of basepairs to both sides',
                    epilog='Uses pybedtools slop and filter')
parser.add_argument('-b', '--bases', action="store", default=50, help = "number of basepairs to add on each side")
parser.add_argument('-i', '--input', action="store", help = ".bed file of repeat locations")
parser.add_argument('-g', '--genome', action="store", help = "genome file (as described in bedtools slop)")
parser.add_argument('-o', '--output', action="store", help = "output file, where to output the slopped .bed file")


args = parser.parse_args()
print("parsing args")
input_features = BedTool(args.input)
print(len(input_features))
expanded_features = input_features.slop(b=int(args.bases), g=args.genome)
expanded_features.saveas(args.output)

    


