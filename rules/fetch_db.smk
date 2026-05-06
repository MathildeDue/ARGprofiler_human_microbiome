rule fetch_db_panres_fa:
    output:
        "prerequisites/db_panres/panres_genes.fa"
    params:
        local_fa=config["panres_fa"],
        time=config["time_path"]
    threads: 1
    log:
        "prerequisites/db_panres/panres_genes.log"
    resources:
        shell_exec="sh"
    shell:
        """
        mkdir -p prerequisites/db_panres

        {params.time} -v --output=prerequisites/db_panres/fetch_panres_fa.bench \
        ln -sf $(realpath {params.local_fa}) {output} > {log} 2>&1
        """


rule fetch_db_panres_meta:
    input:
        fa="prerequisites/db_panres/panres_genes.fa"
    output:
        meta="prerequisites/db_panres/panres_annotations.tsv",
        glengths="prerequisites/db_panres/panres_lengths.tsv"
    params:
        local_meta=config["panres_meta"],
        time=config["time_path"]
    threads: 1
    log:
        "prerequisites/db_panres/panres_meta.log"
    shell:
        """
        mkdir -p prerequisites/db_panres

        {params.time} -v --output=prerequisites/db_panres/fetch_panres_meta.bench \
        ln -sf $(realpath {params.local_meta}) {output.meta} > {log} 2>&1

        awk '
            /^>/ {{
                if (seq_len) print name "\\t" seq_len;
                name = substr($0, 2);
                seq_len = 0;
                next;
            }}
            {{
                seq_len += length($0);
            }}
            END {{
                if (name) print name "\\t" seq_len;
            }}
        ' {input.fa} > {output.glengths}
        """

rule index_db_panres:
    input:
        "prerequisites/db_panres/panres_genes.fa"
    output:
        check_file_index="prerequisites/db_panres/check_file_index_db_panres.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    conda: "../env/kma.yaml"
    params:
        time=config["time_path"]
    threads: 1
    log:
        "prerequisites/db_panres/panres_index.log"
    shell:
        """
        {params.time} -v --output=prerequisites/db_panres/index_panres.bench \
        kma index -i {input} -o prerequisites/db_panres/panres 2> {log}
        touch {output.check_file_index}
        """
        
rule fetch_db_mOTUs:
    output:
        "prerequisites/db_motus/db_mOTU_v3.0.1.tar.gz"
    params:
        zenodo_url=config["zenodo_motus_tar"],
        time=config["time_path"],
    threads: 1
    log:
        "prerequisites/db_motus/db_mOTU_v3.0.1.log"
    shell:
        """
        {params.time} -v --output=prerequisites/db_motus/fetch_mOTUs.bench wget {params.zenodo_url} -P prerequisites/db_motus >> {log}
        """

rule index_db_mOTUs:
    input:
        "prerequisites/db_motus/db_mOTU_v3.0.1.tar.gz"
    output:
        check_file_index="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    conda: 
        "../env/kma.yaml"
    params:
        time=config["time_path"],
    threads: 20
    log:
        "prerequisites/db_motus/index_db_mOTUs.log"
    shell:
        """
        tar -xf {input} -C prerequisites/db_motus/ db_mOTU/db_mOTU_DB_CEN.fasta > {log}
        {params.time} -v --output=prerequisites/db_motus/index_mOTUs.bench kma index -i prerequisites/db_motus/db_mOTU/db_mOTU_DB_CEN.fasta -o prerequisites/db_motus/db_mOTUs 2>> {log}
        touch {output.check_file_index} 
        """
