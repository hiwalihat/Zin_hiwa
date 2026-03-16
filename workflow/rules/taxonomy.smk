# =============================================================================
# Rule module: Taxonomic Classification
# Tools: Kraken2 (k-mer based classification) + Bracken (abundance estimation)
# =============================================================================

rule kraken2_classify:
    """Classify clean reads with Kraken2."""
    input:
        r1 = rules.bowtie2_host_removal.output.r1_clean,
        r2 = rules.bowtie2_host_removal.output.r2_clean,
    output:
        report  = f"{OUT}/taxonomy/{{sample}}/{{sample}}.kraken2.report",
        output  = f"{OUT}/taxonomy/{{sample}}/{{sample}}.kraken2.output",
    params:
        db         = config["databases"]["kraken2_db"],
        confidence = config["kraken2"]["confidence"],
        extra      = config["kraken2"]["extra_args"],
    threads: config["kraken2"]["threads"]
    conda: "../../envs/taxonomy.yaml"
    log: f"{OUT}/logs/taxonomy/{{sample}}_kraken2.log"
    shell:
        """
        kraken2 \
            --db {params.db} \
            --paired {input.r1} {input.r2} \
            --threads {threads} \
            --confidence {params.confidence} \
            --report {output.report} \
            --output {output.output} \
            {params.extra} \
            > {log} 2>&1
        """


rule bracken_abundance:
    """Re-estimate taxonomic abundances from Kraken2 report using Bracken."""
    input:
        report = rules.kraken2_classify.output.report,
    output:
        bracken = f"{OUT}/taxonomy/{{sample}}/{{sample}}.bracken",
        report  = f"{OUT}/taxonomy/{{sample}}/{{sample}}.bracken_report",
    params:
        db        = config["databases"]["bracken_db"],
        read_len  = config["bracken"]["read_length"],
        level     = config["bracken"]["taxonomic_level"],
        threshold = config["bracken"]["threshold"],
    conda: "../../envs/taxonomy.yaml"
    log: f"{OUT}/logs/taxonomy/{{sample}}_bracken.log"
    shell:
        """
        bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} \
            > {log} 2>&1
        """


rule kraken2_mpa_report:
    """Convert Bracken report to MetaPhlAn-style format for downstream use."""
    input:
        report = rules.bracken_abundance.output.report,
    output:
        mpa = f"{OUT}/taxonomy/{{sample}}/{{sample}}.mpa_style.txt",
    conda: "../../envs/taxonomy.yaml"
    log: f"{OUT}/logs/taxonomy/{{sample}}_mpa.log"
    shell:
        """
        kreport2mpa.py \
            -r {input.report} \
            -o {output.mpa} \
            --display-header \
            > {log} 2>&1
        """
