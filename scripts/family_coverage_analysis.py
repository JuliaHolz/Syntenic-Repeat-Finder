import pandas as pd
import subprocess

family_names = ["SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"]

def get_repeat_number(family_name, bed_file):
    cmd = "grep '{};' {} | wc -l".format(family_name, bed_file)
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    return int(ps.communicate()[0])


def quantify_family_through_steps(family_name):
    original_number = get_repeat_number(family_name,"outputs/nosimple_hg38.bed")
    repeats_both_sides = get_repeat_number(family_name,"outputs/f50_hg38.bed") #after_filtering_for_nearby_repeats
    edges_removed = get_repeat_number(family_name,"outputs/6_f50_hg38_e100.bed") #after_expanding_and_removing_too_near_to_edges
    mapped = get_repeat_number(family_name,"outputs/mapped_6_f50_hg38_e100_pantro.bed")
    total_number_mapped = get_repeat_number(family_name, "outputs/pantro_hg38_all_repeats_mapped.bed")
    total_number_unmapped = get_repeat_number(family_name, "outputs/pantro_hg38_all_repeats_unmapped.bed")
    percent_mapped = 0
    if (total_number_mapped+ total_number_unmapped) >0:
        percent_mapped = float(total_number_mapped) / (total_number_mapped+ total_number_unmapped)
    df_list = [family_name, original_number,original_number-repeats_both_sides,repeats_both_sides-edges_removed,edges_removed-mapped, total_number_mapped , total_number_unmapped, percent_mapped]
    return df_list

def build_filter_dataframe(family_list):
    lst = []
    for family in family_list:
        lst.append(quantify_family_through_steps(family))
    df = pd.DataFrame(lst, columns = ["family", "original_number", "repeats_too_close_both_sides", "edge_too_close", "unmapped", "total_mapped", "total_unmapped", "percentage_mapped"])
    return df

def get_family_coverage_dataframe(family_name):
    file_name = "outputs/6_f50_hg38_e100_pantro/alignments/" + family_name + "/repeat_alignment_coverage.csv"
    section_col_names= ["bp", "mis", "indel", "un"]
    left_col_names = ["left"+string for string in section_col_names]
    repeat_col_names = ["repeat"+string for string in section_col_names]
    right_col_names = ["right"+string for string in section_col_names]
    columns = ["id"] + left_col_names + repeat_col_names + right_col_names
    df = pd.read_csv(file_name,names=columns)
    return df

family_list = ["AluYb11", "AluYb8a1", "LTR26C", "LTR26D", "LTR5_Hs", "SVA_E", "SVA_F", "AluYi6_4d","AluYj4","AluYh7"]

more_families = ["ajax", "amalthea", "ananke","aoede", "callisto", "COMP-subunit_5SRNA_rnd-6_family-13719", "COMP-subunit_5SRNA_rnd-6_family-13720", "COMP-subunit_ACRO_rnd-5_family-1624","COMP-subunit_ACRO_rnd-5_family-1625", "COMP-subunit_ACRO_rnd-5_family-37", "COMP-subunit_ACRO_rnd-5_family-38", "COMP-subunit_FAM90A_rnd-6_family-7382", "COMP-subunit_TAF11_rnd-6_family-27360", "COMP-subunit_TELO_rnd-6_family-10479", "COMP-subunit_TELO_rnd-6_family-166", "COMP-subunit_VNTR_rnd-6_family-8746", "COMP-subunit_VNTR_rnd-6_family-8747", "cyllene", "DNM1r", "elara", "erinome", "FAM90Ar", "ghimalia", "harpalyke", "HSAT5v1", "HSAT5v2", "kalyke", "MER5A1r", "pasiphae", "SATR1v", "sinpoe", "SST1v", "teucerv1_5edge", "teucerv2_3edge", "teucerv3_internal", "TIFr", "Walusat", "SN5"]

#df = build_filter_dataframe(family_list+more_families)
#df.to_csv("families_through_filters.csv")


#AluY_shared_families = ["AluY", "AluYa5", "AluYa8", "AluYb8", "AluYb9", "AluYc", "AluYc3", "AluYd8", "AluYe5", "AluYe6", "AluYf1", "AluYg6", "AluYh3", "AluYh9", "AluYi6", "AluYk3", "AluYk4", "AluYk11","AluYk12", "AluYm1"]
#AluY_human_only = ["AluYb11", "AluYb8a1", "AluYi6_4d","AluYj4","AluYh7"]

LTR26_shared_families = ["LTR26", "LTR26B", "LTR26E"]
LTR26_human_only_families = ["LTR26C", "LTR26D"]
df = build_filter_dataframe(LTR26_shared_families+LTR26_human_only_families)
print(df)

print("Shared families")
for family_name in []:

    file_name = "outputs/6_f50_hg38_e100_pantro/alignments/" + family_name + "/repeat_alignment_coverage.csv"
    section_col_names= ["bp", "mis", "indel", "un"]
    left_col_names = ["left"+string for string in section_col_names]
    repeat_col_names = ["repeat"+string for string in section_col_names]
    right_col_names = ["right"+string for string in section_col_names]

    columns = ["id"] + left_col_names + repeat_col_names + right_col_names
    df = pd.read_csv(file_name,names=columns)
    #df = df[df["repeatbp"]>100]
    print(family_name + ": ", df.loc[:, 'repeatun'].mean())
    print(df.head(50))

