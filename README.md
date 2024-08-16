This repository contains code to align repeats (and their surrounding regions) from a target genome to their corresponding regions in a query genome as a way to determine whether the repeat families in the target are shared between the target and query.

To run the full pipeline run the following snakemake command:
snakemake --sdm conda -p outputs/f50_hg38_e100_mLoxAfr_anchorwave/all_alignment_summary.txt
replacing "50" with your filter distance, hg38 with your target name, 100 with the expand distance, mLoxAfr with your query and anchorwave with your aligner (anchorwave or fastga, though fastga is not currently working)

Then, you can go to the notebook in /scripts/detecting_shared_families.ipynb and input the location of the output folder and the clade numbers of your two species, and set the thresholds for considering the families to be shared to explore the number of families the pipeline output matches the repeatmasker annotations for.


The file structure is as follows:
/scripts contains various python scripts (and an ipynb) used by the pipeline
/Snakefile is a snakemake workflow with rules defining how to generate our output files (look here for more documentation on specific steps/rules of the pipeline and software dependencies)
/envs contains a .yml file specifying the conda environment used to run some of the scripts
/inputs contains input files (details below) as well as some files generated from the inputs such as fasta indexes (.fai) and files generated by fastga (.gix and .gdb).
/outputs contains the output files (details below)


Input Files:
To run a new species pair, you need to create the following files:
- /inputs/TARGET_NAME.fna : the full genome of your target species
- /inputs/TARGET_NAME.gff : the full genome annotation (gff file) of your target species
- /inputs/TARGET_NAME.out : the repeatmasker output (.out) when run on your target species
- /inputs/QUERY_NAME.fa : the full genome of your query species


Output Folders:
- /outputs/chains contains chain files and the files used to generate those chain files (psl and sizes.txt files come from using the UCSC pipeline for generating .chain files)
- /outputs/anchorwave alignments contains the whole genome alignments generated by maf in the files QUERY_TARGET.maf, as well as the files generated during the alignment process.
- /outputs/FASTGA_alignments contains the files generated by the FASTGA alignment process
- /outputs/FASTGAtmp contains temporary files for FASTGA
- in the base output directory there will also be all the .bed files generated by the sorting/filtering steps of the pipeline

Output for a specific species pair will be contained in a folder with the following naming convention (items in parentheses will be filled in with specific values/paramenters):
f(filter distance)_(target name)_e(expand distance)_(query name)
Where:
- filter distance is the amount of unique (non-repeat) sequence we require on at least one side of the repeat for the repeat to be included in our analysis
- target name is the name of the target species (must be the same as a .fna, .gff, and .out file in inputs)
- query name is the name of the query species (must be the same as a .fna file in inputs)
- expand distance is the amount of surrounding sequence on each side of the repeat we want to align in base pairs

Within the folder for a specific species pair, there are several files and directories:
Files:
- all_alignment_summary.txt contains the time it took for each family's alignment
- family_summary.txt contains the abundances of each family (how many alignments were actually done, not the original abundances from the repeatmasker output, but the abundances after filtering)
- num_csv_lines_summary.txt contains the total number of lines in the coverage csvs generated for this run (mostly useful for checking that the pipeline ran correctly on all the repeats)
- original_target_intervals_corresponding_to_expanded_and_mapped.bed contains the original intervals from the .out file which are not expanded by the expand distance which make it through the filters and map to the query (due to repeats close to edges you can't always just simply subtract from/add to the end to get the original intervals)

Alignments:
- alignments contains named folders for each family
    - each family folder contains:
        - algnment_summary.txt, which summarizes how long the alignments took
        - alignments.caf, a file containing all the alignments in .caf format, as well as headers with the number of alignments for each repeat
        - errors.txt which has any error output from the aligner
        - repeat_alignment_coverage.csv

repeat_alignment_coverage.csv columns
- id, with format: family name;family classification;repeatmasker id of repeat;line of the original target bed file that the repeat came from
- bases on the left side of the repeat
- percentage of left bases mismatched
- percentage of left bases indeled
- percentage of left bases unaligned
- bases in the repeat
- percentage of repeat bases mismatched
- percentage of repeat bases indeled
- percentage of repeat bases unaligned
- bases on the right side of the repeat
- percentage of right bases mismatched
- percentage of right bases indeled
- percentage of repeat bases unaligned

