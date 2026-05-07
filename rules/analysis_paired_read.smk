rule download_paired_end_reads:
	"""
	Downloading metagenomic raw paired end reads from ENA using enaDataGet
	"""
	output:
		"results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz",
		"results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz",
		check_file_raw="results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt"
	envmodules:
		"tools",
		"fastq-dl/2.0.4",
	conda: "../env/download.yaml"
	params:
		time=config["time_path"],
		attempts=config["max_attempts"]
	threads: 20
	log:
		"results/raw_reads/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/raw_reads/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench fastq-dl -a {wildcards.paired_reads} --silent --cpus {threads} --max-attempts {params.attempts} -o results/raw_reads/paired_end/{wildcards.paired_reads} > {log}
		touch {output.check_file_raw}
		"""

# rule trim_paired_end_reads_fastp:
# 	"""
# 	Adapter trimming of raw paired end reads using fastp, including polyG trimming and deduplication
# 	"""
# 	input:
# 		in1=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz"),
# 		in2=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz")
# 	output:
# 		out1="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
# 		out2="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq",
# 		singleton="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq",
# 		check_file_trim="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt"
# 	params:
# 		overlap_diff_limit=config["overlap_diff_limit"],
# 		average_qual=config["average_qual"],
# 		length_required=config["length_required"],
# 		cut_tail=config["cut_tail"],
# 		poly_g_min_len=config["poly_g_min_len"],
# 		dup_calc_accuracy=config["dup_calc_accuracy"],
# 		out_merge="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_merged.trimmed.fastq",
# 		h="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}.html",
# 		j="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}.json",
# 		time=config["time_path"]
# 	envmodules:
# 		"tools",
# 		"fastp/0.23.2",
# 	conda: "../env/qc.yaml"
# 	threads: 8
# 	log:
# 		"results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}.log"

# 	shell:
# 		"""
# 		mkdir -p results/trimmed_reads_fastp/paired_end/{wildcards.paired_reads}

# 		{params.time} -v --output=results/trimmed_reads_fastp/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
# 		fastp \
# 			-i {input.in1} -I {input.in2} \
# 			-o {output.out1} -O {output.out2} \
# 			--merge --merged_out {params.out_merge} \
# 			--unpaired1 {output.singleton} --unpaired2 {output.singleton} \
# 			--overlap_diff_limit {params.overlap_diff_limit} \
# 			--average_qual {params.average_qual} \
# 			--length_required {params.length_required} \
# 			--trim_poly_g --poly_g_min_len {params.poly_g_min_len} \
# 			--dedup --dup_calc_accuracy {params.dup_calc_accuracy} \
# 			{params.cut_tail} \
# 			-h {params.h} -j {params.j} -w {threads} \
# 			2> {log}

# 		cat {params.out_merge} >> {output.singleton} 2>> {log}
# 		rm {params.out_merge}
# 		touch {output.check_file_trim}
# 		"""

rule trim_paired_end_reads_bbduk:
    """
    Quality and adapter trimming of raw paired end reads using BBDuk
    """
    input:
        in1=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz"),
        in2=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz")
    output:
        out1="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
        out2="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq",
        singleton="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq",
        check_file_trim="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt"
    params:
        ref=config["bbduk_adapter_ref"],
        ktrim=config["bbduk_ktrim"],
        k=config["bbduk_k"],
        mink=config["bbduk_mink"],
        hdist=config["bbduk_hdist"],
        qtrim=config["bbduk_qtrim"],
        trimq=config["bbduk_trimq"],
        minlength=config["bbduk_minlength"],
        trimpolyg=config["bbduk_trimpolyg"],
        overwrite=config["bbduk_overwrite"],
        time=config["time_path"]
    envmodules:
        "openjdk/23.0.1",
        "bbmap/38.90"
    conda:
        "../env/bbduk.yaml"
    threads: 10
    log:
        "results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p results/trimmed_reads_bbduk/paired_end/{wildcards.paired_reads}

        {params.time} -v --output=results/trimmed_reads_bbduk/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        bbduk.sh \
            in1={input.in1} in2={input.in2} \
            out1={output.out1} out2={output.out2} \
            outs={output.singleton} \
            ref={params.ref} \
            ktrim={params.ktrim} k={params.k} mink={params.mink} hdist={params.hdist} \
            qtrim={params.qtrim} trimq={params.trimq} minlength={params.minlength} \
            trimpolyg={params.trimpolyg} \
            tbo tpe threads={threads} overwrite={params.overwrite} \
            2> {log}

        touch {output.check_file_trim}
        """

rule fastp_dup_eval_after_bbduk_before_dedup_paired_end_reads:
    """
    Evaluate duplication on BBDuk-trimmed paired-end reads before deduplication
    """
    input:
        in1="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
        in2="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"
    output:
        out1=temp("results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}_1.eval.fastq"),
        out2=temp("results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}_2.eval.fastq"),
        html="results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}.html",
        json="results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}.json",
        check="results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}_check.txt"
    params:
        dup_calc_accuracy=config["dup_calc_accuracy"],
        time=config["time_path"]
    envmodules:
        "tools",
        "fastp/0.23.2"
    conda:
        "../env/qc.yaml"
    threads: 8
    log:
        "results/fastp_dup_eval/before_bbduk/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p results/fastp_dup_eval/before_bbduk/{wildcards.paired_reads}

        {params.time} -v --output=results/fastp_dup_eval/before_bbduk/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        fastp \
            -i {input.in1} -I {input.in2} \
            -o {output.out1} -O {output.out2} \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            --dup_calc_accuracy {params.dup_calc_accuracy} \
            -h {output.html} -j {output.json} -w {threads} \
            > {log} 2>&1

        touch {output.check}
        """

rule fastp_dup_eval_before_raw_paired_end_reads:
    """
    Evaluate duplication on raw paired-end reads before deduplication
    """
    input:
        in1=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz"),
        in2=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz")
    output:
        out1=temp("results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}_1.eval.fastq"),
        out2=temp("results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}_2.eval.fastq"),
        html="results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}.html",
        json="results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}.json",
        check="results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}_check.txt"
    params:
        dup_calc_accuracy=config["dup_calc_accuracy"],
        time=config["time_path"]
    envmodules:
        "tools",
        "fastp/0.23.2"
    conda:
        "../env/qc.yaml"
    threads: 8
    log:
        "results/fastp_dup_eval/before_raw/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p results/fastp_dup_eval/before_raw/{wildcards.paired_reads}

        {params.time} -v --output=results/fastp_dup_eval/before_raw/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        fastp \
            -i {input.in1} -I {input.in2} \
            -o {output.out1} -O {output.out2} \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            --dup_calc_accuracy {params.dup_calc_accuracy} \
            -h {output.html} -j {output.json} -w {threads} \
            > {log} 2>&1

        touch {output.check}
        """


rule fastp_dedup_bbduk_paired_end_reads:
    """
    Deduplicate BBDuk-trimmed paired-end reads and singleton reads with fastp,
    without further trimming/filtering
    """
    input:
        in1=ancient("results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq"),
        in2=ancient("results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"),
        singleton=ancient("results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq")
    output:
        out1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
        out2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq",
        singleton="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq",
        html="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.html",
        json="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.json",
        check="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check.txt"
    params:
        dup_calc_accuracy=config["dup_calc_accuracy"],
        time=config["time_path"]
    envmodules:
        "tools",
        "fastp/0.23.2"
    conda:
        "../env/qc.yaml"
    threads: 8
    log:
        "results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p results/trimmed_dedup/paired_end/{wildcards.paired_reads}

        # Deduplicate paired-end reads
        {params.time} -v --output=results/trimmed_dedup/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        fastp \
            -i {input.in1} -I {input.in2} \
            -o {output.out1} -O {output.out2} \
            --dedup \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            --dup_calc_accuracy {params.dup_calc_accuracy} \
            -h {output.html} -j {output.json} -w {threads} \
            > {log} 2>&1

        # Deduplicate singleton reads
        fastp \
            -i {input.singleton} \
            -o {output.singleton} \
            --dedup \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            --dup_calc_accuracy {params.dup_calc_accuracy} \
            >> {log} 2>&1

        touch {output.check}
        """

rule fastp_dup_eval_after_dedup_paired_end_reads:
    """
    Evaluate duplication on fastp-deduplicated paired-end reads
    """
    input:
        in1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
        in2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq"
    output:
        out1=temp("results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}_1.eval.fastq"),
        out2=temp("results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}_2.eval.fastq"),
        html="results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}.html",
        json="results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}.json",
        check="results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}_check.txt"
    params:
        dup_calc_accuracy=config["dup_calc_accuracy"],
        time=config["time_path"]
    envmodules:
        "tools",
        "fastp/0.23.2"
    conda:
        "../env/qc.yaml"
    threads: 8
    log:
        "results/fastp_dup_eval/after_dedup/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p results/fastp_dup_eval/after_dedup/{wildcards.paired_reads}

        {params.time} -v --output=results/fastp_dup_eval/after_dedup/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        fastp \
            -i {input.in1} -I {input.in2} \
            -o {output.out1} -O {output.out2} \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            --dup_calc_accuracy {params.dup_calc_accuracy} \
            -h {output.html} -j {output.json} -w {threads} \
            > {log} 2>&1

        touch {output.check}
        """

# rule fastp_dup_eval_after_fastp_before_dedup_paired_end_reads:
#     """
#     Evaluate duplication on fastp trim paired-end reads
#     """
#     input:
#         in1="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
#         in2="results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"
#     output:
#         out1=temp("results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}_1.eval.fastq"),
#         out2=temp("results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}_2.eval.fastq"),
#         html="results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}.html",
#         json="results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}.json",
#         check="results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}_check.txt"
#     params:
#         dup_calc_accuracy=config["dup_calc_accuracy"],
#         time=config["time_path"]
#     envmodules:
#         "tools",
#         "fastp/0.23.2"
#     conda:
#         "../env/qc.yaml"
#     threads: 8
#     log:
#         "results/fastp_dup_eval/before_fastp/{paired_reads}/{paired_reads}.log"
#     shell:
#         """
#         mkdir -p results/fastp_dup_eval/before_fastp/{wildcards.paired_reads}

#         {params.time} -v --output=results/fastp_dup_eval/before_fastp/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
#         fastp \
#             -i {input.in1} -I {input.in2} \
#             -o {output.out1} -O {output.out2} \
#             --disable_adapter_trimming \
#             --disable_quality_filtering \
#             --disable_length_filtering \
#             --disable_trim_poly_g \
#             --dup_calc_accuracy {params.dup_calc_accuracy} \
#             -h {output.html} -j {output.json} -w {threads} \
#             > {log} 2>&1

#         touch {output.check}
#         """


rule fastqc_raw_paired_end_reads:
	"""
	Run FastQC on raw paired-end reads
	"""
	input:
		read_1=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz"),
		read_2=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz")
	output:
		html_1="results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}_1_fastqc.html",
		zip_1="results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}_1_fastqc.zip",
		html_2="results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}_2_fastqc.html",
		zip_2="results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}_2_fastqc.zip",
		check_file_fastqc="results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}_check_file_fastqc.txt"
	envmodules:
		"tools",
		"fastqc"
	conda: "../env/qc.yaml"
	params:
		outdir="results/fastqc/raw/paired_end/{paired_reads}",
		time=config["time_path"]
	threads: 4
	log:
		"results/fastqc/raw/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		mkdir -p {params.outdir}
		{params.time} -v --output=results/fastqc/raw/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
		fastqc -t {threads} -o {params.outdir} {input.read_1} {input.read_2} > {log} 2>&1
		touch {output.check_file_fastqc}
		"""

# rule fastqc_fastp_paired_end_reads:
# 	"""
# 	Run FastQC on fastp-trimmed paired-end reads
# 	"""
# 	input:
# 		read_1=ancient("results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq"),
# 		read_2=ancient("results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"),
# 		read_3=ancient("results/trimmed_reads_fastp/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq")
# 	output:
# 		html_1="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_1.trimmed_fastqc.html",
# 		zip_1="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_1.trimmed_fastqc.zip",
# 		html_2="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_2.trimmed_fastqc.html",
# 		zip_2="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_2.trimmed_fastqc.zip",
# 		html_3="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed_fastqc.html",
# 		zip_3="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed_fastqc.zip",
# 		check_file_fastqc="results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}_check_file_fastqc.txt"
# 	envmodules:
# 		"tools",
# 		"fastqc"
# 	conda: "../env/qc.yaml"
# 	params:
# 		outdir="results/fastqc/fastp/paired_end/{paired_reads}",
# 		time=config["time_path"]
# 	threads: 4
# 	log:
# 		"results/fastqc/fastp/paired_end/{paired_reads}/{paired_reads}.log"
# 	shell:
# 		"""
# 		mkdir -p {params.outdir}
# 		{params.time} -v --output=results/fastqc/fastp/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
# 		fastqc -t {threads} -o {params.outdir} {input.read_1} {input.read_2} {input.read_3} > {log} 2>&1
# 		touch {output.check_file_fastqc}
# 		"""

rule fastqc_dedup_bbduk_paired_end_reads:
    """
    Run FastQC on BBDuk-trimmed + fastp-deduplicated paired-end reads
    """
    input:
        read_1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
        read_2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq",
        read_3="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq"
    output:
        html_1="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup_fastqc.html",
        zip_1="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup_fastqc.zip",
        html_2="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup_fastqc.html",
        zip_2="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup_fastqc.zip",
        html_3="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup_fastqc.html",
        zip_3="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup_fastqc.zip",
        check_file_fastqc="results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check_file_fastqc.txt"
    envmodules:
        "tools",
        "fastqc"
    conda:
        "../env/qc.yaml"
    params:
        outdir="results/fastqc/trimmed_dedup/paired_end/{paired_reads}",
        time=config["time_path"]
    threads: 4
    log:
        "results/fastqc/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p {params.outdir}
        {params.time} -v --output=results/fastqc/trimmed_dedup/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        fastqc -t {threads} -o {params.outdir} {input.read_1} {input.read_2} {input.read_3} > {log} 2>&1
        touch {output.check_file_fastqc}
        """

rule kma_paired_end_reads_mOTUs:
    """
    Mapping deduplicated paired reads for identifying AMR using KMA with mOTUs db
    """
    input:
        read_1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
        read_2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq",
        read_3="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq",
        check_file_db_mOTUs="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
    output:
        "results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.res",
        "results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.mapstat",
        check_file_kma_mOTUs="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt"
    params:
        db="prerequisites/db_motus/db_mOTUs",
        outdir="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}",
        kma_params=config["kma_params_motus"],
        time=config["time_path"]
    envmodules:
        "tools",
        "kma/1.4.12a",
    conda:
        "../env/kma.yaml"
    threads: 20
    log:
        "results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.log"
    shell:
        """
        {params.time} -v --output=results/kma_mOTUs/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        kma -ipe {input.read_1} {input.read_2} -i {input.read_3} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2>> {log}
        rm -f results/kma_mOTUs/paired_end/{wildcards.paired_reads}/*.aln 2>> {log} || true

        if [ -f results/kma_mOTUs/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.fsa ]; then
            gzip -f results/kma_mOTUs/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.fsa 2>> {log}
        fi
        touch {output.check_file_kma_mOTUs}
        """

rule kma_paired_end_reads_panRes:
	"""
	Mapping raw paired reads for identifying AMR using KMA with panres db
	"""
	input: 
		read_1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
		read_2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq",
		read_3="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq",
		check_file_db_panres="prerequisites/db_panres/check_file_index_db_panres.txt"
	output:
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.res",
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mat.gz",
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat",
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bam",
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered",
		check_file_kma_panres="results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt"
	params:
		db="prerequisites/db_panres/panres",
		outdir="results/kma_panres/paired_end/{paired_reads}/{paired_reads}",
		kma_params=config["kma_params_panres"],
		mapstat="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat",
		mapstat_filtered="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered",
		mapstat_table="prerequisites/db_panres/panres_lengths.tsv",
		time=config["time_path"]
	envmodules:
		"tools",
		"kma/1.4.12a",
		"samtools/1.16",
		"gcc/9.4.0",
		"intel/perflibs/64/2020_update2",
		"R/4.3.0",
	threads: 2
	conda: "../env/kma_panres.yaml"
	log:
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/kma_panres/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench kma -ipe {input.read_1} {input.read_2} -i {input.read_3} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2> {log} |samtools fixmate -m - -|samtools view -u -bh -F 4|samtools sort -o results/kma_panres/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bam 2>> {log}
		rm -f results/kma_panres/paired_end/{wildcards.paired_reads}/*.aln 2>> {log} || true

		if [ -f results/kma_panres/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.fsa ]; then
			gzip -f results/kma_panres/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.fsa 2>> {log}
		fi

		Rscript prerequisites/mapstat_filtering/mapstatFilters.R \
			-i {params.mapstat} \
			-o {params.mapstat_filtered} \
			-r {params.mapstat_table} 2>> {log}	
		touch {output.check_file_kma_panres}
		"""

rule mash_sketch_paired_end_reads:
	"""
	Creation of mash sketches of paired end reads using mash
	"""
	input:
		read_1=ancient("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq"),
		read_2=ancient("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq"),
		read_3=ancient("results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq"),
	output:
		out="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.dedup.fastq.msh",
		check_file_mash="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt"
	envmodules:
		"tools",
		"mash/2.3",
	conda: "../env/mash.yaml"
	params:
		time=config["time_path"],
		k=config["mash_k"],
		s=config["mash_s"]
	threads: 20
	log:
		"results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/mash_sketch/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench cat {input.read_1} {input.read_2} {input.read_3} | mash sketch -k {params.k} -s {params.s} -I {wildcards.paired_reads} -C Paired -r -o {output.out} -p {threads} - 2>> {log}
		touch {output.check_file_mash}
		"""

rule mash_triangle_paired_end_reads:
	"""
	Compare all paired-end mash sketches
	"""
	input:
		lambda wildcards: expand(
			"results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.dedup.fastq.msh",
			paired_reads=paired
		)
	output:
		out="results/mash_sketch/paired_end/mash_triangle.tsv"
	params:
		time=config["time_path"]
	envmodules:
		"tools",
		"mash/2.3"
	conda: "../env/mash.yaml"
	threads: 20
	log:
		"results/mash_sketch/paired_end/mash_triangle.log"
	shell:
		"""
		{params.time} -v --output=results/mash_sketch/paired_end/mash_triangle.bench \
		mash triangle {input} > {output.out} 2> {log}
		"""

rule mash_dist_paired_end_reads:
	"""
	Create a combined mash sketch file and compare all paired-end sketches all-vs-all with mash dist
	"""
	input:
		lambda wildcards: expand(
			"results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.dedup.fastq.msh",
			paired_reads=paired
		)
	output:
		mash_all="results/mash_sketch/paired_end/all_samples.msh",
		mash_tab="results/mash_sketch/paired_end/mash_dist.tsv",
		matrix="results/mash_sketch/paired_end/mash_distance_matrix.tsv"
	params:
		time=config["time_path"],
		script="scripts/create_distance_matrix.sh",
		mash_prefix="results/mash_sketch/paired_end/all_samples"
	envmodules:
		"tools",
		"mash/2.3"
	conda: "../env/mash.yaml"
	threads: 20
	log:
		"results/mash_sketch/paired_end/mash_dist.log"
	shell:
		"""
		set -euo pipefail

		{params.time} -v --output=results/mash_sketch/paired_end/mash_paste.bench \
		mash paste {params.mash_prefix} {input} 2> {log}

		{params.time} -v --output=results/mash_sketch/paired_end/mash_dist.bench \
		mash dist {output.mash_all} {output.mash_all} > {output.mash_tab} 2>> {log}

		bash {params.script} {output.mash_tab} {output.matrix} 2>> {log}
		"""

rule sourmash_sketch_dedup_paired_end_reads:
	"""
	Create per-sample sourmash signatures from paired end deduplicated reads
	"""
	input:
		read_1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
		read_2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq",
		read_3="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq"
	output:
		out="results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.sig",
		check_file_sourmash="results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt"
	params:
		tmp1="results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_r1.sig",
		tmp2="results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_r2.sig",
		tmp3="results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.sig",
		scaled=1000,
		k=31,
		time=config["time_path"]
	envmodules:
		"tools",
		"sourmash"
	conda: "../env/sourmash.yaml"
	threads: 20
	log:
		"results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		mkdir -p results/sourmash/trimmed_dedup/paired_end/{wildcards.paired_reads}
		{params.time} -v --output=results/sourmash/trimmed_dedup/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
		sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp1} {input.read_1} > {log} 2>&1
		sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp2} {input.read_2} >> {log} 2>&1
		sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp3} {input.read_3} >> {log} 2>&1
		sourmash sig merge {params.tmp1} {params.tmp2} {params.tmp3} -o {output.out} --name {wildcards.paired_reads} >> {log} 2>&1
		rm -f {params.tmp1} {params.tmp2} {params.tmp3}
		touch {output.check_file_sourmash}
		"""

rule sourmash_compare_dedup_paired_end_reads:
	"""
	Compare all paired-end sourmash signatures
	"""
	input:
		lambda wildcards: expand(
			"results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}.sig",
			paired_reads=paired
		)
	output:
		matrix="results/sourmash/trimmed_dedup/paired_end/sourmash_matrix.npy",
		csv="results/sourmash/trimmed_dedup/paired_end/sourmash_distances.csv"
	params:
		time=config["time_path"]
	envmodules:
		"tools",
		"sourmash"
	conda: "../env/sourmash.yaml"
	threads: 20
	log:
		"results/sourmash/trimmed_dedup/paired_end/sourmash_compare.log"
	shell:
		"""
		{params.time} -v --output=results/sourmash/trimmed_dedup/paired_end/sourmash_compare.bench \
		sourmash compare {input} -o {output.matrix} --csv {output.csv} > {log} 2>&1
		"""

rule sourmash_sketch_trimmed_paired_end_reads:
    """
    Create per-sample sourmash signatures from BBDuk-trimmed paired-end reads
    """
    input:
        read_1="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
        read_2="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq",
        read_3="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq"
    output:
        out="results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}.sig",
        check_file_sourmash="results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt"
    params:
        tmp1="results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}_r1.sig",
        tmp2="results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}_r2.sig",
        tmp3="results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}_singleton.sig",
        scaled=1000,
        k=31,
        time=config["time_path"]
    envmodules:
        "tools",
        "sourmash"
    conda: "../env/sourmash.yaml"
    threads: 20
    log:
        "results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}.log"
    shell:
        """
        mkdir -p results/sourmash/trimmed_bbduk/paired_end/{wildcards.paired_reads}

        {params.time} -v --output=results/sourmash/trimmed_bbduk/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
        sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp1} {input.read_1} > {log} 2>&1

        sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp2} {input.read_2} >> {log} 2>&1

        sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp3} {input.read_3} >> {log} 2>&1

        sourmash sig merge {params.tmp1} {params.tmp2} {params.tmp3} -o {output.out} --name {wildcards.paired_reads} >> {log} 2>&1

        rm -f {params.tmp1} {params.tmp2} {params.tmp3}
        touch {output.check_file_sourmash}
        """

rule sourmash_compare_trimmed_paired_end_reads:
    """
    Compare all paired-end trimmed sourmash signatures
    """
    input:
        lambda wildcards: expand(
            "results/sourmash/trimmed_bbduk/paired_end/{paired_reads}/{paired_reads}.sig",
            paired_reads=paired
        )
    output:
        matrix="results/sourmash/trimmed_bbduk/paired_end/sourmash_matrix.npy",
        csv="results/sourmash/trimmed_bbduk/paired_end/sourmash_distances.csv"
    params:
        time=config["time_path"]
    envmodules:
        "tools",
        "sourmash"
    conda: "../env/sourmash.yaml"
    threads: 20
    log:
        "results/sourmash/trimmed_bbduk/paired_end/sourmash_compare.log"
    shell:
        """
        {params.time} -v --output=results/sourmash/trimmed_bbduk/paired_end/sourmash_compare.bench \
        sourmash compare {input} -o {output.matrix} --csv {output.csv} > {log} 2>&1
        """
    
rule sourmash_sketch_raw_paired_end_reads:
	"""
	Create per-sample sourmash signatures from raw paired end reads
	"""
	input:
		read_1=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz"),
		read_2=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz")
	output:
		out="results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}.sig",
		check_file_sourmash="results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt"
	params:
		tmp1="results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}_r1.sig",
		tmp2="results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}_r2.sig",
		scaled=1000,
		k=31,
		time=config["time_path"]
	envmodules:
		"tools",
		"sourmash"
	conda: "../env/sourmash.yaml"
	threads: 20
	log:
		"results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		mkdir -p results/sourmash/raw/paired_end/{wildcards.paired_reads}

		{params.time} -v --output=results/sourmash/raw/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
		sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp1} {input.read_1} > {log} 2>&1

		sourmash sketch dna -p scaled={params.scaled},k={params.k},abund -o {params.tmp2} {input.read_2} >> {log} 2>&1

		sourmash sig merge {params.tmp1} {params.tmp2} -o {output.out} --name {wildcards.paired_reads} >> {log} 2>&1

		rm -f {params.tmp1} {params.tmp2}
		touch {output.check_file_sourmash}
		"""

rule sourmash_compare_raw_paired_end_reads:
	"""
	Compare all paired-end raw sourmash signatures
	"""
	input:
		lambda wildcards: expand(
			"results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}.sig",
			paired_reads=paired
		)
	output:
		matrix="results/sourmash/raw/paired_end/sourmash_matrix.npy",
		csv="results/sourmash/raw/paired_end/sourmash_distances.csv"
	params:
		time=config["time_path"]
	envmodules:
		"tools",
		"sourmash"
	conda: "../env/sourmash.yaml"
	threads: 20
	log:
		"results/sourmash/raw/paired_end/sourmash_compare.log"
	shell:
		"""
		{params.time} -v --output=results/sourmash/raw/paired_end/sourmash_compare.bench \
		sourmash compare {input} -o {output.matrix} --csv {output.csv} > {log} 2>&1
		"""

rule ARG_extender_paired_reads:
    """
    Performing local ARG extension of paired reads using perl script
    """
    input:
        read_1="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_1.dedup.fastq",
        read_2="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_2.dedup.fastq",
        read_3="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_singleton.dedup.fastq",
        panres_mapstat_filtered="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered"
    output:
        out_fasta="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.fasta.gz",
        out_gfa="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.gfa.gz",
        out_frag="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag.gz",
        out_frag_gz="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag_raw.gz",
        check_file_ARG="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}_check_file_ARG.txt"
    params:
        # The number of iterations for the ARG extender. "-1" means that the extension will continue until no new contigs are found.
        ARG="25",
        temp_dir="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}",
        db="prerequisites/db_panres/panres_genes.fa",
        time=config["time_path"]
    conda: "../env/argextender.yaml"
    threads: 20
    log:
        "results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.log"
    shell:
        """
        if grep -q -v -m 1 "#" {input.panres_mapstat_filtered};
        then
            echo "running argextender" > {log}
            {params.time} -v --output=results/ARG_extender/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench \
            perl prerequisites/ARGextender/targetAsm.pl {params.ARG} {threads} {params.temp_dir} {params.db} {input.read_1} {input.read_2} {input.read_3} 2>> {log}
            gzip -f results/ARG_extender/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.fasta 2>> {log}
            gzip -f results/ARG_extender/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.gfa 2>> {log}
            touch {output.check_file_ARG}
        else
            echo "not running argextender" > {log}
            touch {output.out_fasta}
            touch {output.out_gfa}
            touch {output.out_frag}
            touch {output.out_frag_gz}
            touch {output.check_file_ARG}
        fi
        """


rule cleanup_paired_end_reads:
    """
    Removing unwanted files
    """
    input:
        check_file_raw="results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt",
        check_file_trim="results/trimmed_reads_bbduk/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt",
        check_file_dedup="results/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check.txt",
        check_file_mash="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt",
        check_file_sourmash_trim="results/sourmash/trimmed_dedup/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt",
        check_file_sourmash_raw="results/sourmash/raw/paired_end/{paired_reads}/{paired_reads}_check_file_sourmash.txt",
        check_file_kma_mOTUs="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt",
        check_file_kma_panres="results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt",
        check_file_ARG="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}_check_file_ARG.txt"
    output:
        check_file_clean_final1="results/raw_reads/paired_end/{paired_reads}/check_clean_raw.txt",
        check_file_clean_final2="results/trimmed_reads_bbduk/paired_end/{paired_reads}/check_clean_trim.txt"
    shell:
        """
        rm -f results/raw_reads/paired_end/{wildcards.paired_reads}/*.gz
        touch {output.check_file_clean_final1}
        touch {output.check_file_clean_final2}
        """