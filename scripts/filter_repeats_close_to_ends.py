import argparse 
from pybedtools import BedTool

parser = argparse.ArgumentParser(
                    prog='FilterRepeatsCloseToEnds',
                    description='takes in a sorted bed file of genomic locations, and filters out repeats close to the ends, now no longer used in pipeline',
                    epilog='Uses pybedtools  filter')
parser.add_argument('-b', '--bases', action="store", default=50, help = "number of basepairs added on each side with slop")
parser.add_argument('-i', '--input', action="store", help = "bed file containing 6 columns for the slopped/expanded repeats with 6 columns for the non-expanded repeats to the right")
parser.add_argument('-o', '--output', action="store", help = "location to output the bed file to")


args = parser.parse_args()


input_features = BedTool(args.input)
expected_diff = 2*int(args.bases)
#check that the difference in sizes is what we expect -- 
filtered_features = input_features.filter(lambda i:((int(i.fields[2])-int(i.fields[1])) - (int(i.fields[8])-int(i.fields[7]))==expected_diff))

filtered_features.saveas(args.output)
