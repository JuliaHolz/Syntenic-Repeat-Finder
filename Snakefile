import os
wildcard_constraints:
    genome = "hg38",
    filterdist = "[\\d]+",
    bp = "[\\d]+",
    target = "[^_]*",
    query = "[^_]*"

rule repeatmasker_output_to_bed:
    input: 
        repeatmasker_out = "inputs/{target}.out",
        script = "scripts/RMout_to_bed.py"
    output:  "outputs/{target}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -i {input.repeatmasker_out} -o {output}"

#First we remove the simple/low complexity regions
rule filter_simple_repeats_and_low_complexity:
    input: "outputs/{filename}.bed"
    output: "outputs/nosimple_{filename}.bed"
    shell: "awk -F'[;\\t]'  '{{ if(($5 !~ /Simple_repeat.*/) &&($5 !~ /Low_complexity.*/)) print }}' {input} > {output}"


rule sort_bed:
    input: "outputs/nosimple_{filepath}.bed"
    output: "outputs/{filepath}_sorted.bed"
    conda: "envs/pybedtools.yml"
    shell: "bedtools sort -i {input} > {output}"

rule filter_for_adjacent:
    input: 
        repeats = "outputs/{target}_sorted.bed",
        script = "scripts/filter_for_adjacent_features.py"
    output: "outputs/f{filterdist}_{target}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -d {wildcards.filterdist} -i {input.repeats} -f {input.repeats} -o {output}"

rule index_genome:
    input: "inputs/{target}.fna"
    output: "outputs/{target}.fna.fai"
    conda: "envs/pybedtools.yml"
    shell: "/u3/local/samtools/bin/samtools faidx {input} -o {output}"

rule create_genome_file_for_bedtools:
    input: "outputs/{target}.fna.fai"
    output: "outputs/{target}_genomeFile.txt"
    conda: "envs/pybedtools.yml"
    shell: "awk -v OFS='\\t' {{'print $1,$2'}} {input} > {output}"

#uses bedtools slop to add basepairs on both sides, not going past the ends of the chromosomes
rule add_basepairs_on_both_sides:
    input:
        repeats = "outputs/f{filterdist}_{target}.bed",
        genome = "outputs/{target}_genomeFile.txt",
        script = "scripts/add_basepairs_on_each_side.py"
    output: "outputs/f{filterdist}_{target}_e{bp}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -b {wildcards.bp} -i {input.repeats} -g {input.genome} -o {output}"

#combine the bed files side by side so we can filter out the ones that were too close to the start/end to slop outwards
#rule combine_beds:
#    input:
#        not_slopped = "outputs/f{filterdist}_{target}.bed",
#        slopped = "outputs/f{filterdist}_{target}_slop_{bp}.bed"
#    output: "outputs/f{filterdist}_{target}_combined_{bp}.bed"
#    conda: "envs/pybedtools.yml"
#    shell: "paste -d\"\\t\" {input.slopped} {input.not_slopped}>{output}"

#rule filter_repeats_close_to_ends:
#    input:
#        repeats =  "outputs/f{filterdist}_{target}_combined_{bp}.bed",
#        script = "scripts/filter_repeats_close_to_ends.py"
#    output: "outputs/f{filterdist}_{target}_e{bp}.bed"
#    conda: "envs/pybedtools.yml"
#    shell: "python {input.script} -b {wildcards.bp} -i {input.repeats} -o {output}"

#rule remove_last_6_cols_of_bed:
#    input: "outputs/f{filterdist}_{target}_e{bp}.bed"
#    output: "outputs/6_f{filterdist}_{target}_e{bp}.bed"
#    conda: "envs/pybedtools.yml"
#    shell: "cut -d$'\t' -f 1-6 {input} > {output}"

rule lift_repeats:
    input:
        chain = "inputs/{target}_{query}.chain",
        repeats = "outputs/f{filterdist}_{target}_e{bp}.bed"
    output: 
        mapped = "outputs/mapped_f{filterdist}_{target}_e{bp}_{query}.bed",
        unmapped = "outputs/unmapped_f{filterdist}_{target}_e{bp}_{query}.bed"
    conda: "envs/pybedtools.yml"
    shell: "/usr/local/ucscTools/liftOver {input.repeats} {input.chain} {output.mapped} {output.unmapped}"

#gets a bed file of only the repeats that mapped over in human
rule get_corresponding_unmapped_repeats:
    input:
        mapped_repeats = "outputs/mapped_f{filterdist}_{target}_e{bp}_{query}.bed",
        target_repeats = "outputs/f{filterdist}_{target}_e{bp}.bed"
    output: "outputs/orig_corresponding_to_mapped_f{filterdist}_{target}_e{bp}_{query}.bed"
    shell: "awk -F'\\t' 'NR==FNR{{c[$4]++;next}};c[$4] > 0' {input.mapped_repeats} {input.target_repeats} > {output}"




#split the files into beds by family, note, we use gensub to replace any / in family names with % so we can use them as file names
checkpoint split_file_by_families:
    input: 
        mapped_repeats = "outputs/mapped_{repeatfile}.bed",
        corresponding_repeats = "outputs/orig_corresponding_to_mapped_{repeatfile}.bed"
    output: "outputs/{repeatfile}/family_summary.txt"
    run:
        shell("mkdir -p outputs/{wildcards.repeatfile}")
        shell("mkdir -p outputs/{wildcards.repeatfile}/query_beds")
        shell("mkdir -p outputs/{wildcards.repeatfile}/target_beds")
        shell("awk -F'[;\\t]' '{{print>(\"outputs/{wildcards.repeatfile}/query_beds/\" gensub(\"/\", \"%\", \"g\", $4) \".bed\")}}' {input.mapped_repeats}")
        shell("awk -F'[;\\t]' '{{print>(\"outputs/{wildcards.repeatfile}/target_beds/\" gensub(\"/\", \"%\", \"g\", $4) \".bed\")}}' {input.corresponding_repeats}")
        shell("awk -F '[\\t;]' '{{print $4}}' {input.mapped_repeats} | sort | uniq -c > {output}")
        shell("mkdir -p outputs/{wildcards.repeatfile}/alignments")
        shell("mkdir -p outputs/{wildcards.repeatfile}/target_fasta")
        shell("mkdir -p outputs/{wildcards.repeatfile}/query_fasta")



rule align_family:
    input: 
        mapped_bed = "outputs/f{filterdist}_{target}_e{bp}_{query}/query_beds/{family}.bed", 
        target_bed = "outputs/f{filterdist}_{target}_e{bp}_{query}/target_beds/{family}.bed",
        target_genome = "inputs/{target}.fna",
        query_genome = "inputs/{query}.fna", 
        script = "scripts/align_family.py"
    output: 
        query_fasta = "outputs/f{filterdist}_{target}_e{bp}_{query}/query_fasta/{family}.fasta",
        target_fasta = "outputs/f{filterdist}_{target}_e{bp}_{query}/target_fasta/{family}.fasta",
        summary = "outputs/f{filterdist}_{target}_e{bp}_{query}/alignments/{family}/alignment_summary.txt",
        caf = "outputs/f{filterdist}_{target}_e{bp}_{query}/alignments/{family}/alignments.caf"
    threads: 1
    conda: "envs/pybedtools.yml"
    shell: """python {input.script} -i {input.target_bed} -m {input.mapped_bed} -t {input.target_genome} -q {input.query_genome} -o outputs/f{wildcards.filterdist}_{wildcards.target}_e{wildcards.bp}_{wildcards.query} -f {wildcards.family}"""


def aggregate_families(wildcards):
     checkpoint_output = checkpoints.split_file_by_families.get(**wildcards).output[0]
     sections = checkpoint_output.split("/")
     path = sections[0] + "/" + sections[1] + "/query_beds/"
     families, = glob_wildcards(os.path.join(path, "{family}.bed"))
     path_to_summary = sections[0] + "/" + sections[1] + "/alignments/"
     return expand(os.path.join(path_to_summary, "{FAMILY}/alignment_summary.txt"), FAMILY=families)


rule align_all_families:
    input: 
        aggregate_families
    output: "outputs/{repeatfile}/all_alignment_summary.txt"
    shell: "cat ./outputs/{wildcards.repeatfile}/alignments/*/alignment_summary.txt > {output}"

checkpoint target_interval_bed_per_family:
    input:
        target_bed = "outputs/f{filterdist}_{target}.bed",
        query_bed = "outputs/mapped_f{filterdist}_{target}_e{bp}_{query}.bed"
    output: 
        all_target = "outputs/f{filterdist}_{target}_e{bp}_{query}/original_target_intervals_corresponding_to_expanded_and_mapped.bed",
    run:
        shell("mkdir -p outputs/f{wildcards.filterdist}_{wildcards.target}_e{wildcards.bp}_{wildcards.query}/nonexpanded_target_beds")
        shell("awk -F'\\t' 'NR==FNR{{c[$4]++;next}};c[$4] > 0' {input.query_bed} {input.target_bed} > {output.all_target}")
        shell("awk -F'[;\\t]' '{{print>(\"outputs/f{wildcards.filterdist}_{wildcards.target}_e{wildcards.bp}_{wildcards.query}/nonexpanded_target_beds/\" gensub(\"/\", \"%\", \"g\", $4) \".bed\")}}' {output.all_target}")




#remove the 7_19_ once I rerun the rest of the pipeline
rule parse_family_caf:
    input:
        script = "scripts/pared_down_parse_family_caf.py",
        target_beds_done = "outputs/f{filterdist}_{target}_e{bp}_{query}/original_target_intervals_corresponding_to_expanded_and_mapped.bed",
        family_target_bed = "outputs/f{filterdist}_{target}_e{bp}_{query}/nonexpanded_target_beds/{family}.bed",
        family_expanded_bed = "outputs/f{filterdist}_{target}_e{bp}_{query}/target_beds/{family}.bed",
        family_caf = "outputs/f{filterdist}_{target}_e{bp}_{query}/alignments/{family}/alignments.caf"
    output:
        family_table = "outputs/f{filterdist}_{target}_e{bp}_{query}/alignments/{family}/repeat_alignment_coverage.csv"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -c {input.family_caf} -r {input.family_target_bed} -e {input.family_expanded_bed} -o {output}"

def aggregate_tsvs(wildcards):
     checkpoint_output = checkpoints.target_interval_bed_per_family.get(**wildcards).output[0]
     sections = checkpoint_output.split("/")
     path = sections[0] + "/" + sections[1] + "/query_beds/"
     families, = glob_wildcards(os.path.join(path, "{family}.bed"))
     path_to_summary = sections[0] + "/" + sections[1] + "/alignments/"
     return expand(os.path.join(path_to_summary, "{FAMILY}/repeat_alignment_coverage.csv"), FAMILY=families)


rule parse_all_caf: 
    input: 
        aggregate_tsvs
    output: "outputs/f{filterdist}_{target}_e{bp}_{query}/num_csv_lines_summary.txt"
    threads: 1
    shell: "cat ./outputs/f{wildcards.filterdist}_{wildcards.target}_e{wildcards.bp}_{wildcards.query}/alignments/*/repeat_alignment_coverage.csv | wc -l > {output}"
