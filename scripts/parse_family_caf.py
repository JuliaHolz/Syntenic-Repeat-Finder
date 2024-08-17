import subprocess
import argparse 
import time
from pybedtools import BedTool
import re
from enum import Enum



parser = argparse.ArgumentParser(
                    prog='ParseFamilyCaf',
                    description='takes in a file consisting of a series of .caf alignments output by crossmatch for a repeat family and parses out the information we are interested in'
                )

parser.add_argument('-c', '--caf', action="store", help = "file containing a series of .caf alignments")
parser.add_argument('-r', '--repeats', action = "store", help = "bed file containing locations of the repeats in the target genome that are found in the alignments")
parser.add_argument('-e', '--expanded', action = "store", help = "bed file containing locations of the repeats (after being expanded to contain flanking regions on each side) in the target genome that are found in the alignments")
parser.add_argument('-o', '--output', action = "store", help = "location for outputting the table")

args = parser.parse_args()


class Column(Enum):
    match = 1
    mismatch = 2
    deleted = 3
    unaligned = 4

def columnToChar(col:Column):
    if col == Column.match:
        return '.'
    if col == Column.mismatch:
        return 'M'
    if col == Column.deleted:
        return '-'
    if col == Column.unaligned:
        return '~'
    raise Exception("invalid column " + str(col))

def isACTG(char):
    return (char != "-" and char != "+" and char != "/" and char != ",")
def get_section_stat_string(coverage):
    percent_mismatch = (float(coverage.count(Column.mismatch))/len(coverage)) * 100
    percent_deleted = (float(coverage.count(Column.deleted))/len(coverage)) * 100 
    percent_unaligned = (float(coverage.count(Column.unaligned))/len(coverage)) * 100
    return "{},{:.2f},{:.2f},{:.2f}".format(len(coverage),percent_mismatch,percent_deleted,percent_unaligned)

def get_csv_line(coverage, repeat_bedline, expanded_bedline):
    left_flanking_to_cover = repeat_bedline.start - expanded_bedline.start
    repeat_end_idx = repeat_bedline.end - expanded_bedline.start
    left_flanking_coverage = coverage[0:left_flanking_to_cover]
    repeat_coverage = coverage[left_flanking_to_cover:repeat_end_idx]
    right_flanking_coverage = coverage[repeat_end_idx:]
    id = repeat_bedline.name.replace(",", "#")
    return id+","+ get_section_stat_string(left_flanking_coverage) + "," + get_section_stat_string(repeat_coverage) + "," + get_section_stat_string(right_flanking_coverage) + "\n"



#process alignments computes coverage for the repeat and surrounding area
#if a base is matched in any alignment it is considered matched, the next priority is mismatch, then deleted, 
#and a base will only be considered unaligned if it is unaligned in all alignments
def processAlignments(repeat_id, alignment_list, repeat_bedline, expanded_bedline):
    left_flanking_to_cover = repeat_bedline.start - expanded_bedline.start
    right_flanking_to_cover = expanded_bedline.end - repeat_bedline.end
    coverage = [Column.unaligned for i in range(0,expanded_bedline.end-expanded_bedline.start)]
    for alignment in alignment_list:
        #check that this is an alignment of the repeat to itself in the other genome
        alignment_list = alignment.split(",")
        qid = alignment_list[4].split(":")[1]
        #check that this is an alignment of the repeat we're looking for -- if this fails some alignments or ids are not being recognized correctly
        assert(qid.strip()==repeat_id.replace("/", "%").strip())
        #subtract 1 because the caf indices are fully closed, one indexed
        qstart = int(alignment_list[5])
        coverage_idx = qstart-1
        alignment_string = alignment_list[16]
        in_qdeleted = False
        in_insertion = False
        for char_idx in range(0, len(alignment_string)):
            #first -- normal case -- letter is a match 
            if in_qdeleted:
                if(alignment_string[char_idx] == '-'): #end the indel
                    in_qdeleted = False 
                elif(isACTG(alignment_string[char_idx])):
                    #this is a nucleotide that is deleted in the query -- only change if unaligned - mismatch and match are higher "ranking" 
                    if(coverage[coverage_idx] == Column.unaligned):
                        coverage[coverage_idx] = Column.deleted
                    coverage_idx +=1
                else:
                    raise Exception("Invalid character in deletion section: " + str(alignment_string[char_idx]))
            elif in_insertion:
                if(alignment_string[char_idx] == '+'): #end the insertion
                    in_insertion = False 
                elif not isACTG(alignment_string[char_idx]):
                    raise Exception("Invalid character in insertion section: " + str(alignment_string[char_idx]))
            elif alignment_string[char_idx] == '+':
                in_insertion = True
            elif alignment_string[char_idx] == '-':
                in_qdeleted = True
            elif alignment_string[char_idx] == '/':
                #mismatch should be two non-special chars
                assert(char_idx>0 and char_idx+1<len(alignment_string) and isACTG(alignment_string[char_idx-1]) and isACTG(alignment_string[char_idx+1]))
            elif (char_idx+1 < len(alignment_string)) and alignment_string[char_idx+1] == '/':
                if(coverage[coverage_idx] != Column.match): #anything other than match is lower "rank" than mismatch
                    coverage[coverage_idx] = Column.mismatch
                coverage_idx +=1
            elif (char_idx> 0) and alignment_string[char_idx-1] == '/':
                pass #we already dealt with mismatch when we encountered the first char of it
            elif isACTG(alignment_string[char_idx]): #match case
                coverage[coverage_idx] = Column.match
                coverage_idx +=1
            else:
                raise Exception("Invalid character:" + str(alignment_string[char_idx]) + " or malformed caf alignment string") 
        qend = int(alignment_list[6])
        assert(coverage_idx == qend)
    #debugging - print out string representing coverage of the expanded interval
    #try:
    #    coverage_string = "".join([columnToChar(cov) for cov in coverage])
    #    print(coverage_string)
    #except:
    #    print(coverage)
    assert(len(coverage)==expanded_bedline.end-expanded_bedline.start)
    return get_csv_line(coverage, repeat_bedline, expanded_bedline)

bed_idx = 0
caf_id = ""
id_found = False
repeat_id = ""
num_alignments_found = False
num_alignments_expected = -1
num_alignments_found = -1


target_repeat = BedTool(args.repeats)
expanded_target_repeat = BedTool(args.expanded)
alignment_list = []
num_alignments_expected = 0
with open(args.caf, "r") as caf_file,open(args.output, "w") as outfile:
    for line in caf_file:
        if line[0] == ">":
            #id lines look like this: >1 alignments for:5S#rRNA#29137#34877
            repeat_id = line.split(":")[1]
            repeat_bedline = target_repeat[bed_idx]
            split_bedline = repeat_bedline.name.split(";")
            expanded_bedline = expanded_target_repeat[bed_idx]
            bed_idx+=1
            #check that IDs are the same -- input bedfiles should have the same order of entries as the alignments.caf file
            repeat_id_tuple= (split_bedline[2], split_bedline[3])
            expanded_id_tuple =  (expanded_bedline.name.split(";")[2], expanded_bedline.name.split(";")[3])
            alignment_id_tuple = (repeat_id.split("#")[2], repeat_id.split("#")[3].strip())
            assert(repeat_id_tuple==expanded_id_tuple)
            assert(repeat_id_tuple==alignment_id_tuple)
            #determine how many alignments we expect
            num_alignments_expected = int(line.split()[0][1:])
            alignment_list = []
            if num_alignments_expected == 0:
                outfile.write(processAlignments(repeat_id, alignment_list, repeat_bedline,expanded_bedline))
        else:
            alignment_list.append(line)
            if (len(alignment_list)==num_alignments_expected):
                outfile.write(processAlignments(repeat_id, alignment_list, repeat_bedline,expanded_bedline))
            


