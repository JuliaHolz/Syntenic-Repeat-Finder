import re
import argparse 


parser = argparse.ArgumentParser(
                    prog='RMoutToBed',
                    description='takes in a repeatmasker output file, and creates a bed file',
                    epilog='Uses re to match repeatmasker output lines with regular expressions')
parser.add_argument('-i', '--input', action="store", help = "repeatmasker .out file")
parser.add_argument('-o', '--output', action="store", help = "location to output the bed file to")
args = parser.parse_args()

repeatmasker_output_file = open(args.input, 'r')
bed_file = open(args.output, "w")

#initial regex taken from robert's script
#location_line_expression = re.compile(r"^\s*(?P<score>\d+)\s+(?P<perc_div>\d+.\d)+\s+\d+.\d+\s+\d+.\d+.*")
#r means raw string so we don't have to escape the backslashes
full_location_line_expression = re.compile(r"^\s*(?P<score>\d+)\s+(?P<perc_div>\d+.\d+)\s+(?P<perc_del>\d+.\d+)\s+(?P<perc_ins>\d+.\d+)\s+(?P<query_sequence>\S+)\s+(?P<query_position_begin>\d+)\s+(?P<query_position_end>\d+)\s+(?P<query_left>\(\d+\))\s+(?P<orientation>[+C])\s+(?P<matching_repeat>\S+)\s+(?P<class_family>\S+)\s+(?P<repeat_begin>\(?-?\d+\)?)\s+(?P<repeat_end>\(?-?\d+\)?)\s+(?P<repeat_left>\(?\d+\)?)\s+(?P<id>\d+).*")


count = 0
print_parsed_attributes = False
line_number = -1
while True:
 
    # Get next line from file
    line = repeatmasker_output_file.readline()
    line_number += 1
    # if line is empty
    # end of file is reached
    if not line:
        break
    #if we match the expected format of the line we can use this line to make a bed output file
    matched = full_location_line_expression.match(line)
    if matched != None:
        if(print_parsed_attributes):
            print("MATCHED LINE")
            print("score:", matched.group("score"))
            print("perc div:", matched.group("perc_div"))
            print("perc del:", matched.group("perc_del"))
            print("perc ins:", matched.group("perc_ins"))
            print("query sequence:", matched.group("query_sequence"))
            print("query position begin:", matched.group("query_position_begin"))
            print("query position end:", matched.group("query_position_end"))
            print("query left:", matched.group("query_left"))
            print("orientation:", matched.group("orientation"))
            print("matched repeat:", matched.group("matching_repeat"))
            print("class family:", matched.group("class_family"))

            print("repeat position begin:", matched.group("repeat_begin"))
            print("repeat left:", matched.group("repeat_left"))
            print("ID:", matched.group("id"))
        orientation = matched.group("orientation")
        strandString = "-" if orientation=="C" else ("+" if orientation=="+" else ".")
        #output string has standard columns of bed6 format, the name column is the name of the matching repeat, the class/family, the repeatmasker id, and the line number it's found on in the bed file (zero-indexed)
        name_string = matched.group("matching_repeat")+";"+ matched.group("class_family")+";"+matched.group("id")+";"+str(line_number)
        outputString = matched.group("query_sequence") + "\t" +  matched.group("query_position_begin") + "\t" + matched.group("query_position_end") +"\t" + name_string +"\t" + matched.group("score")  + "\t" + strandString + "\n"
        bed_file.write(outputString)
        count+=1
            
        if count % 100000 == 0:
            print(str(count) + " lines written to file")
        
    else:
        print("Unmatched line in file, printed below, should be a non-element line:")
        print(line)
        print()
            

        

repeatmasker_output_file.close()
bed_file.close()