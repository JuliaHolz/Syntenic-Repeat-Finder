import argparse 
from pybedtools import BedTool
#This script takes in a  bed file of repeat locations, and filters out those without unique sequence on one or more sides, outputting a new, filtered bed file

#options when running this script
#-d, minimum distance from other features (what is considered "too close") in basepairs
#-i, input file -- this should be the .bed file of repeat locations, which need to all have distinct names for filtering to work appropriately (N=True option)
#-f, feature file -- this is the set of features we want to check the distance to (may, but does not have to be the same as -i file)
#-o, output file -- where to output the filtered .bed file
#-t, by default, the filter is one sided, requiring at least d distance from the nearest repeat on one side, but if you use the -t flag it will require d distance on both sides
parser = argparse.ArgumentParser(
                    prog='FilterAdjacent',
                    description='takes in a sorted bed file of genomic locations, and filters out those without unique sequence on one or more sides, outputting a new, filtered bed file.',
                    epilog='Uses pybedtools closest and filter')
parser.add_argument('-d', '--distance', action="store", default=50, help = "minimum distance from other repeats (what is considered \"too close\") in basepairs")
parser.add_argument('-i', '--input', action="store", help = ".bed file of repeat locations, which need to all have distinct names for filtering to work appropriately")
parser.add_argument('-f', '--features', action="store", help = ".bed file of repeats we want to check the distance to (may, but does not have to be the same as -i file)")
parser.add_argument('-o', '--output', action="store", help = "output file, where to output the filtered .bed file")
parser.add_argument('-t', '--twosided', action="store_const",const=True, help = "by default, the filter is one sided, requiring at least d distance from the nearest repeat on one side, but if you use the -t flag it will require d distance on both sides")



args = parser.parse_args()



input_features = BedTool(args.input)
other_features = BedTool(args.features)
sided_tag = ""
d = int(args.distance)
if(args.twosided):
    sided_tag = "t"
    #want to use the -N option for bedtools closest so repeat is not counted as intersecting itself (requires different name for closest feature)
    closest_repeats = input_features.closest(other_features, d=True,t="first", N=True)
    #the last field of the interval is the distance to the nearest other repeat, so we filter on the nearest other repeat being a certain distance away
    closest_repeats_filtered = closest_repeats.filter(lambda i:int(i.fields[-1])>=d)
else:
    #id is ignore downstream, and iu is ignore upstream options
    closest_repeats_upstream = input_features.closest(other_features, D="ref",t="first", N=True, id=True)
    closest_repeats_both_directions = closest_repeats_upstream.closest(other_features, D="ref",t="first", N=True, iu=True)
    #NOTE!!!! MAKE SURE THAT THE INDICES OF THE FIELDS BEING CHECKED ARE CORRECT (BASED ON COLUMNS IN  THE ORIGINAL .BED) AND YOU'RE FILTERING BASED ON DISTANCE, NOT SOME OTHER COLUMN OF THE BED FILE
    #-1 and -8 are based on what fields are in the bed file output by my script for converting repeatmasker output
    closest_repeats_filtered = closest_repeats_both_directions.filter(lambda i:(abs(int(i.fields[-1]))>=d) or (abs(int(i.fields[-8]))>=d))

closest_repeats_filtered.bed6().saveas(args.output)