# =============================================================================
# Rule module: Host Read Removal
# Tool: Bowtie2  – align to host genome, keep unmapped (non-host) reads
# =============================================================================

rule bowtie2_host_removal:
    """Align trimmed reads against host genome; discard mapped (host) reads."""
    input:
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2,
    output:
        r1_clean = f"{OUT}/host_removed/{{sample}}/{{sample}}_R1_clean.fastq.gz",
        r2_clean = f"{OUT}/host_removed/{{sample}}/{{sample}}_R2_clean.fastq.gz",
        stats    = f"{OUT}/host_removed/{{sample}}/{{sample}}_bowtie2.log",
    params:
        index  = config["databases"]["host_genome"],
        extra  = config["bowtie2"]["extra_args"],
    threads: config["bowtie2"]["threads"]
    conda: "../../envs/host_removal.yaml"
    log: f"{OUT}/logs/host_removal/{{sample}}.log"
    shell:
        """
        bowtie2 \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            {params.extra} \
            --un-conc-gz {OUT}/host_removed/{wildcards.sample}/{wildcards.sample}_clean_%.fastq.gz \
            -S /dev/null \
            2> {output.stats}

        # Rename outputs to expected names
        mv {OUT}/host_removed/{wildcards.sample}/{wildcards.sample}_clean_1.fastq.gz \
           {output.r1_clean}
        mv {OUT}/host_removed/{wildcards.sample}/{wildcards.sample}_clean_2.fastq.gz \
           {output.r2_clean}

        cp {output.stats} {log}
        """
