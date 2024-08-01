import os

#To run this pipeline you need several things
#software:
#conda (to create the pybedtools environment from /envs/pybedtools.yml (which installs pybedtools))
#ucsc tools
#minimap
#LAST (from which we use the maf-convert script)
#AnchorWave 

#inputs:
#target fna: entire target genome, should be located in inputs/target.fna, where target is the name you want in file names relating to this genome (target name may not contain underscores)
#target gff: .gff annotation file for target genome, should be located in inputs/target.gff
#query fna: entire query genome, should be located in inputs/query.fna, where query is the name you want in file names relating to this genome (query name may not contain underscores)
#repeatmasker out file: .out file generated by running repeatmasker on the target genome, should be located in inputs/target.out

#Note: the term "target" here refers to the genome we want to use as the source for our repeats, and "query" refers to the genome we want to lift the repeats over to. Essentially, this pipeline aims to answer this question: which of the repeat families identified by repeatmasker in the target are present at syntenic locations in the query?

#prevent our wildcards from matching _s so snakemake's greedy matching doesn't give us grief
wildcard_constraints:
    genome = "[^_]*",
    filterdist = "[\\d]+",
    bp = "[\\d]+",
    target = "[^_]*",
    query = "[^_]*",

#AnchorWave alignment steps are based on instructions found here: https://github.com/baoxingsong/AnchorWave/blob/master/doc/guideline.pdf
#(the AnchorWave github)
#the anchorwave version I used is anchorwave v1.2.3 with an update to deal with overlapping cds sequences (courtesy of Robert Hubley)
anchorwave_location = "/home/rhubley/projects/multi-species-alignment/AnchorWave-1.2.3/anchorwave"
#the minimap version I used is minimap 2.28-r1209
minimap_location = "/home/jholz/minimap2/minimap2-2.28_x64-linux/minimap2"
#location of LAST for maf-convert
maf_convert_location = "/usr/local/last-1410/bin/maf-convert"
#location of UCSC liftover
ucscTools_location = "/usr/local/ucscTools"

#Generate CDS sequences from .gff file -- this takes some time
rule generate_cds:
    input:
        fna = "inputs/{target}.fna",
        gff = "inputs/{target}.gff"
    output: "outputs/anchorwave_alignments/{target}_cds.fa"
    shell: anchorwave_location + " gff2seq -r {input.fna} -i {input.gff} -o {output}"

#use minimap to align the cds against a genome (we need to do this for both target and query)
rule align_cds_to_genome:
    input:
        genome = "inputs/{genome}.fna",
        cds = "outputs/anchorwave_alignments/{cds_genome}_cds.fa"
    output: "outputs/anchorwave_alignments/{genome}_cds{cds_genome}.sam"
    shell: minimap_location + " -x splice -t 11 -k 12 -a -p 0.4 -N 20 {input.genome} {input.cds} > {output}"

#create the anchors
rule create_anchors:
    input: 
        target_gff = "inputs/{target}.gff",
        target_fna = "inputs/{target}.fna",
        target_cds_alignment = "outputs/anchorwave_alignments/{target}_cds{target}.sam",
        query_cds_alignment = "outputs/anchorwave_alignments/{query}_cds{target}.sam",
        cds = "outputs/anchorwave_alignments/{target}_cds.fa",
        query_fna = "inputs/{query}.fna",
    output: "outputs/anchorwave_alignments/{target}_{query}.anchors"
    shell:  anchorwave_location + " proali -i {input.target_gff} -r {input.target_fna} -a {input.query_cds_alignment} -as {input.cds} -ar {input.query_cds_alignment} -s {input.query_fna} -n {output} -R 1 -Q 1 -ns"    

#do the alignment -- note f_maf is the alignment file for interanchor regions, and maf is the alignment file for the whole genome
rule anchorwave_alignment:
    input:
        target_gff = "inputs/{target}.gff",
        target_fna = "inputs/{target}.fna",
        target_cds_alignment = "outputs/anchorwave_alignments/{target}_cds{target}.sam",
        query_cds_alignment = "outputs/anchorwave_alignments/{query}_cds{target}.sam",
        cds = "outputs/anchorwave_alignments/{target}_cds.fa",
        query_fna = "inputs/{query}.fna",
        anchors = "outputs/anchorwave_alignments/{target}_{query}.anchors"
    output: 
        maf = "outputs/anchorwave_alignments/{target}_{query}.maf",
        f_maf = "outputs/anchorwave_alignments/{target}_{query}.f.maf",
        log =  "outputs/anchorwave_alignments/{target}_{query}.log"
    shell: "/usr/bin/time " + anchorwave_location +" proali -i {input.target_gff} -r {input.target_fna}  -a {input.query_cds_alignment} -as {input.cds} -ar {input.target_cds_alignment} -s {input.query_fna} -n {input.anchors} -o {output.maf} -t 1 -R 1 -Q 1 -B -2 -O1 -4 -E1 -2 -O2 -80 -E2 -1 -f {output.f_maf} -w 38000 -fa3 200000 > {output.log} 2>&1"

#generate the chain file 
#we use LAST's maf-convert script here, though you could also set it up to use UCSC's pipeline (axtChain, chinPreNet, netChainSubset)
#the chain file generation step is a bit of a black box, but it seemed that using UCSC's tools instead of maf-convert led to slightly fewer repeats mapping over, 
#presumably because there is a filtering step in the UCSC pipeline which discards some lower-quality chains

rule generate_chain:
    input: "outputs/anchorwave_alignments/{target}_{query}.maf"
    output: "outputs/chains/{target}_{query}.chain"
    shell: maf_convert_location + " chain {input} > {output}"

#uncomment these rules and comment the generate_chain rule to switch to using the UCSC pipeline for generating the chain file
'''
rule get_chromosome_sizes:
    input: "inputs/{genome}.fna"
    output: "outputs/chains/{genome}_sizes.txt"
    shell: ucscTools_location + "/faSize -detailed {input} > {output}"

rule convert_to_psl:
    input: "outputs/anchorwave_alignments/{target}_{query}.maf"
    output: "outputs/chains/{target}_{query}.psl"
    shell: maf_convert_location + " psl {input} > {output}"

rule axt_chain:
    input: 
        psl = "outputs/chains/{target}_{query}.psl",
        target_fna = "inputs/{target}.fna",
        query_fna = "inputs/{query}.fna"
    output: "outputs/chains/{target}_{query}.axtchain"
    shell: ucscTools_location + "/axtChain -linearGap=loose -psl {input.psl} -faQ -faT {inputs.target_fna} {inputs.query_fna} {output}"

rule chain_prenet: 
    input: 
        target_sizes = "outputs/chains/{target}_sizes.txt",
        query_sizes = "outputs/chains/{query}_sizes.txt",
        axt_chain = "outputs/chains/{target}_{query}.axtchain"
    output: "outputs/chains/{target}_{query}.preNet"
    shell: ucscTools_location + "/chainPreNet {input.axt_chain} {input.target_sizes} {input.query_sizes} {output}"

rule chain_net:
    input: 
        target_sizes = "outputs/chains/{target}_sizes.txt",
        query_sizes = "outputs/chains/{query}_sizes.txt",
        prenet = "outputs/chains/{target}_{query}.preNet"
    output:
        ref_target = "outputs/chains/{target}_{query}_refTarget.chainNet",
        species_target = "outputs/chains/{target}_{query}_speciesTarget.chainNet"
    shell: ucscTools_location + "/chainNet {input.prenet} {input.target_sizes} {input.query_sizes} {output.ref_target} {output.species_target}"

rule final_chain:
    input:
        ref_target =  "outputs/chains/{target}_{query}_refTarget.chainNet",
        prenet = "outputs/chains/{target}_{query}.preNet"
    output: "outputs/chains/{target}_{query}.chain"
    shell: ucscTools_location + "/netChainSubset {input.ref_target} {input.prenet} {output}"
'''

#parses repeatmasker .out file to create a bed file
rule repeatmasker_output_to_bed:
    input: 
        repeatmasker_out = "inputs/{target}.out",
        script = "scripts/RMout_to_bed.py"
    output:  "outputs/{target}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -i {input.repeatmasker_out} -o {output}"

#first we remove the simple/low complexity regions
rule filter_simple_repeats_and_low_complexity:
    input: "outputs/{filename}.bed"
    output: "outputs/nosimple_{filename}.bed"
    shell: "awk -F'[;\\t]'  '{{ if(($5 !~ /Simple_repeat.*/) &&($5 !~ /Low_complexity.*/)) print }}' {input} > {output}"

#sort beds (many bed commands require us to sort our beds)
rule sort_bed:
    input: "outputs/nosimple_{filepath}.bed"
    output: "outputs/{filepath}_sorted.bed"
    conda: "envs/pybedtools.yml"
    shell: "bedtools sort -i {input} > {output}"

#filter out repeats that have other repeats too close to them on both sides
rule filter_for_adjacent:
    input: 
        repeats = "outputs/{target}_sorted.bed",
        script = "scripts/filter_for_adjacent_features.py"
    output: "outputs/f{filterdist}_{target}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -d {wildcards.filterdist} -i {input.repeats} -f {input.repeats} -o {output}"

#we need the index for bedtools to make a genome file
rule index_genome:
    input: "inputs/{target}.fna"
    output: "outputs/{target}.fna.fai"
    conda: "envs/pybedtools.yml"
    shell: "/u3/local/samtools/bin/samtools faidx {input} -o {output}"

#creates a file that is needed by bedtools to perform the "slop" operation (names and lengths of chromosomes) using the samtools index of the genome .fna
rule create_genome_file_for_bedtools:
    input: "outputs/{target}.fna.fai"
    output: "outputs/{target}_genomeFile.txt"
    conda: "envs/pybedtools.yml"
    shell: "awk -v OFS='\\t' {{'print $1,$2'}} {input} > {output}"

#uses bedtools slop to add bp basepairs on both sides of the repeat regions, not going past the ends of the chromosomes
rule add_basepairs_on_both_sides:
    input:
        repeats = "outputs/f{filterdist}_{target}.bed",
        genome = "outputs/{target}_genomeFile.txt",
        script = "scripts/add_basepairs_on_each_side.py"
    output: "outputs/f{filterdist}_{target}_e{bp}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -b {wildcards.bp} -i {input.repeats} -g {input.genome} -o {output}"

#removed these steps because we should be able to deal with repeats that are close to the edges in our coverage calculations
#these steps combine the bed files side by side so we can filter out the ones that were too close to the start/end to slop outwards
'''
rule combine_beds:
    input:
        not_slopped = "outputs/f{filterdist}_{target}.bed",
        slopped = "outputs/f{filterdist}_{target}_slop_{bp}.bed"
    output: "outputs/f{filterdist}_{target}_combined_{bp}.bed"
    conda: "envs/pybedtools.yml"
    shell: "paste -d\"\\t\" {input.slopped} {input.not_slopped}>{output}"

rule filter_repeats_close_to_ends:
    input:
        repeats =  "outputs/f{filterdist}_{target}_combined_{bp}.bed",
        script = "scripts/filter_repeats_close_to_ends.py"
    output: "outputs/f{filterdist}_{target}_e{bp}.bed"
    conda: "envs/pybedtools.yml"
    shell: "python {input.script} -b {wildcards.bp} -i {input.repeats} -o {output}"

rule remove_last_6_cols_of_bed:
    input: "outputs/f{filterdist}_{target}_e{bp}.bed"
    output: "outputs/6_f{filterdist}_{target}_e{bp}.bed"
    conda: "envs/pybedtools.yml"
    shell: "cut -d$'\t' -f 1-6 {input} > {output}"
'''

#lift the repeats from the target genome to the query (requires a .chain file at the location target_query.chain in inputs
rule lift_repeats:
    input:
        chain = "inputs/{target}_{query}.chain",
        repeats = "outputs/f{filterdist}_{target}_e{bp}.bed"
    output: 
        mapped = "outputs/mapped_f{filterdist}_{target}_e{bp}_{query}.bed",
        unmapped = "outputs/unmapped_f{filterdist}_{target}_e{bp}_{query}.bed"
    conda: "envs/pybedtools.yml"
    shell: ucscTools_location + "/liftOver {input.repeats} {input.chain} {output.mapped} {output.unmapped}"

#gets a bed file of only the repeats that mapped over into query in the target (used by our parsing scripts)
rule get_corresponding_unmapped_repeats:
    input:
        mapped_repeats = "outputs/mapped_f{filterdist}_{target}_e{bp}_{query}.bed",
        target_repeats = "outputs/f{filterdist}_{target}_e{bp}.bed"
    output: "outputs/orig_corresponding_to_mapped_f{filterdist}_{target}_e{bp}_{query}.bed"
    shell: "awk -F'\\t' 'NR==FNR{{c[$4]++;next}};c[$4] > 0' {input.mapped_repeats} {input.target_repeats} > {output}"

#split the files into beds by family, note, we use gensub to replace any / in family names with % so we can use them as file names
#this is a checkpoint because we don't know exactly what/how many families we will find in the bed ahead of time
#this script has outputs of a bunch of repeat files for each family -- we don't include these in the outputs because:
# a) there are often 1000+ of them, so snakemake doing modification checking on them makes things slow
# b) we don't know how many there will be (we could give snakemake an output directory but then we run into the issue mentioned in a)
# so please don't mess with the query_beds/family.bed and target_beds/family.bed files or delete the alignments, target_fasta or query_fasta folders after running this step
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


#uses crossmatch (since it is singlethreaded which allows us to run many families at once) to align all instances of a family
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

#need to deal with parens in family names here 
def aggregate_families(wildcards):
     #asking for checkpoint output forces the checkpoint job split_file_by_families to run before this
     checkpoint_output = checkpoints.split_file_by_families.get(**wildcards).output[0]
     sections = checkpoint_output.split("/")
     path = sections[0] + "/" + sections[1] + "/query_beds/"
     families, = glob_wildcards(os.path.join(path, "{family}.bed"))
     path_to_summary = sections[0] + "/" + sections[1] + "/alignments/"
     return expand(os.path.join(path_to_summary, "{FAMILY}/alignment_summary.txt"), FAMILY=families)

#rule to find all the family files, which we don't know how many there will be (hence the function-as-input)
rule align_all_families:
    input: 
        aggregate_families
    output: "outputs/{repeatfile}/all_alignment_summary.txt"
    shell: "cat ./outputs/{wildcards.repeatfile}/alignments/*/alignment_summary.txt > {output}"

# splits the lifted/mapped beds by family for use by our .caf parser
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

#parses in our file of all the alignments (lines with repeat name and alignments in CAF/yaCAF file format)
#outputs a csv with the coverage information for each family
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

#rule to force all the familys' .caf files to be parsed in
rule parse_all_caf: 
    input: 
        aggregate_tsvs
    output: "outputs/f{filterdist}_{target}_e{bp}_{query}/num_csv_lines_summary.txt"
    threads: 1
    shell: "cat ./outputs/f{wildcards.filterdist}_{wildcards.target}_e{wildcards.bp}_{wildcards.query}/alignments/*/repeat_alignment_coverage.csv | wc -l > {output}"

#awk command: awk -F';' 'previd==$3{count1 = split(prevline, prev,  " ", seps); count2 = split($0, curr, " ", seps2); prevline = curr[1] "\t" prev[2] "\t" (prev[3]<curr[3]?curr[3]:prev[3]) "\t" prev[4] "\t" prev[5] "\t" prev[6]} previd!=$3{print prevline; prevline = $0}{previd=$3}' sorted.txt