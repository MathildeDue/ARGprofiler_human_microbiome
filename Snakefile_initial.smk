import sys
import json
import glob
import os
import re
import subprocess
from collections import defaultdict

configfile: "config/config.yaml"


def build_inputs(jsonFile, local_folder="local_reads"):
    
    # Sets to store sample IDs in
    single, paired = set(), set()

    # Open the json file
    with open(jsonFile, "r") as jf:
        data = json.load(jf)

    # Loop through the json input and sort the file ids
    for sample_id, sample_id_read_type in data.items():
        if sample_id_read_type["type"] == "PAIRED":
            paired.add(sample_id)
        elif sample_id_read_type["type"] == "SINGLE":
            single.add(sample_id)

    # Now extract files in local_folder
    local_reads = sorted(glob.glob(os.path.join(local_folder, "*.fastq.gz")))
    single_files, paired_files = defaultdict(set), defaultdict(set)

    # Handles both:
    # sample_1.fastq.gz / sample_2.fastq.gz
    # sample_R1.fastq.gz / sample_R2.fastq.gz
    p = re.compile(r'(.+?)(?:[_]?R|_)([12])\.fastq\.gz$')

    # Now sort the files whether they are paired reads or not
    for fastqFile in local_reads:
        fname = os.path.basename(fastqFile)
        m = p.match(fname)

        if m:
            sample_id = m.group(1)
            paired_files[sample_id].add(fastqFile)
            paired.add(sample_id)
        else:
            sample_id = fname.split(".fastq")[0]
            single_files[sample_id].add(fastqFile)
            single.add(sample_id)

    cmd_paired_template = "ln -sf {} {} && ln -sf {} {} && touch {}_check_file_raw.txt"
   
    #.format(os.path.realpath(fastqFile), os.path.join(dest_dir, os.path.basename(fastqFile)))
    for paired_id, files in paired_files.items():
        # double check there are exactly two files stored for the id
        files = sorted(list(files))
        if len(files) == 2:
            dest_dir = os.path.join("results", "raw_reads", "paired_end", paired_id)
            os.makedirs(dest_dir, exist_ok=True)

            # decide which file is read 1 and read 2
            f1, f2 = files
            b1, b2 = os.path.basename(f1), os.path.basename(f2)

            if "_R1.fastq.gz" in b1 or "_1.fastq.gz" in b1:
                read1, read2 = f1, f2
            else:
                read1, read2 = f2, f1

            cmd = cmd_paired_template.format(
                os.path.realpath(read1),
                os.path.join(dest_dir, f"{paired_id}_1.fastq.gz"),
                os.path.realpath(read2),
                os.path.join(dest_dir, f"{paired_id}_2.fastq.gz"),
                os.path.join(dest_dir, paired_id),
            )
            subprocess.run(cmd, shell=True, check=True)

    cmd_single_template = "ln -sf {} {} && touch {}_check_file_raw.txt"
    
    # move single end files
    for single_id, files in single_files.items():
        files = list(files)

        if len(files) == 1:
            dest_dir = os.path.join("results", "raw_reads", "single_end", single_id)
            os.makedirs(dest_dir, exist_ok=True)

            cmd = cmd_single_template.format(
                os.path.realpath(files[0]),
                os.path.join(dest_dir, os.path.basename(files[0])),
                os.path.join(dest_dir, single_id),
            )
            subprocess.run(cmd, shell=True)

    return list(single), list(paired)



single, paired = build_inputs("input.json")

include: "rules/analysis_paired_read.smk"
include: "rules/analysis_single_read.smk"
include: "rules/fetch_db.smk"


rule all:
    input:
        # Databases
        "prerequisites/db_panres/check_file_index_db_panres.txt",
        "prerequisites/db_panres/panres_lengths.tsv",

        # Raw / trim cleanup checks
        expand("results/raw_reads/paired_end/{paired_reads}/check_clean_raw.txt", paired_reads=paired),
        expand("results/trimmed_reads_bbduk/paired_end/{paired_reads}/check_clean_trim.txt", paired_reads=paired),
        expand("results/raw_reads/single_end/{single_reads}/check_clean_raw.txt", single_reads=single),
        expand("results/trimmed_reads/single_end/{single_reads}/check_clean_trim.txt", single_reads=single),

        # Fastp-trimmed paired-end reads
        expand("results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt", paired_reads=paired),
        expand("results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_check_file_fastqc.txt", paired_reads=paired),

        # BBDuk-trimmed reads
        expand("results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt", paired_reads=paired),
        expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}_check_file_trim.txt", single_reads=single),

        # FastQC
        expand("results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}_check_file_fastqc.txt", paired_reads=paired),
        expand("results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check_file_fastqc.txt", paired_reads=paired),

        # fastp duplication evaluation BEFORE dedup
        expand("results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}_check.txt", paired_reads=paired),
        expand("results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}_check.txt", paired_reads=paired),
        expand("results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}_check.txt", paired_reads=paired),

        # fastp deduplicated outputs
        expand("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq", paired_reads=paired),
        expand("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq", paired_reads=paired),
        expand("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq", paired_reads=paired),
        expand("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check.txt", paired_reads=paired),

        # fastp duplication evaluation AFTER dedup
        expand("results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}_check.txt", paired_reads=paired),

        # KMA mOTUs
        expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
        expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
        expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),

        expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res", single_reads=single),
        expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
        expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),

        # KMA panres
        expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
        expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mat.gz", paired_reads=paired),
        expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
        expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered", paired_reads=paired),
        expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bam", paired_reads=paired),
        expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),

        expand("results/kma_panres/single_end/{single_reads}/{single_reads}.res", single_reads=single),
        expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mat.gz", single_reads=single),
        expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
        expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered", single_reads=single),
        expand("results/kma_panres/single_end/{single_reads}/{single_reads}.bam", single_reads=single),
        expand("results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),

        # Mash
        expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.dedup.fastq.msh", paired_reads=paired),
        expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt", paired_reads=paired),
        expand("results/mash_sketch/single_end/{single_reads}/{single_reads}.trimmed.fastq.msh", single_reads=single),
        expand("results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt", single_reads=single),

        "results/mash_sketch/paired_end/mash_triangle.tsv",
        "results/mash_sketch/paired_end/mash_dist.tsv",
        "results/mash_sketch/paired_end/mash_distance_matrix.tsv",

        # Sourmash trimmed + deduplicated paired-end
        expand("results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.sig", paired_reads=paired),
        expand("results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt", paired_reads=paired),
        "results/sourmash/trimmed_dedup/paired_end/sourmash_matrix.npy",
        "results/sourmash/trimmed_dedup/paired_end/sourmash_distances.csv",

        # Sourmash raw paired-end
        expand("results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}.sig", paired_reads=paired),
        expand("results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt", paired_reads=paired),
        "results/sourmash/raw/paired_end/sourmash_matrix.npy",
        "results/sourmash/raw/paired_end/sourmash_distances.csv",

        # Sourmash BBDuk-trimmed paired-end
        expand("results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}.sig", paired_reads=paired),
        expand("results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt", paired_reads=paired),
        "results/sourmash/trimmed_bbduk/paired_end/sourmash_matrix.npy",
        "results/sourmash/trimmed_bbduk/paired_end/sourmash_distances.csv",

        # ARG extender
        expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.fasta.gz", paired_reads=paired),
        expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.gfa.gz", paired_reads=paired),
        expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag.gz", paired_reads=paired),
        expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag_raw.gz", paired_reads=paired),
        expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}_check_file_ARG.txt", paired_reads=paired),

        expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.fasta.gz", single_reads=single),
        expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.gfa.gz", single_reads=single),
        expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.frag.gz", single_reads=single),
        expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.frag_raw.gz", single_reads=single),
        expand("results/ARG_extender/single_end/{single_reads}/{single_reads}_check_file_ARG.txt", single_reads=single)