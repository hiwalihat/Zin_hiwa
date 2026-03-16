# =============================================================================
# Rule module: Quality Control
# Tools: FastQC (pre-trim), fastp (adapter trimming), FastQC (post-trim)
# =============================================================================

rule fastqc_raw:
    """Run FastQC on raw reads before trimming."""
    input:
        r1 = lambda wc: f"{config['reads_dir']}/{wc.sample}_R1.fastq.gz",
        r2 = lambda wc: f"{config['reads_dir']}/{wc.sample}_R2.fastq.gz",
    output:
        html_r1 = f"{OUT}/qc/fastqc/{{sample}}_R1_fastqc.html",
        html_r2 = f"{OUT}/qc/fastqc/{{sample}}_R2_fastqc.html",
        zip_r1  = f"{OUT}/qc/fastqc/{{sample}}_R1_fastqc.zip",
        zip_r2  = f"{OUT}/qc/fastqc/{{sample}}_R2_fastqc.zip",
    params:
        outdir = f"{OUT}/qc/fastqc",
    threads: 2
    conda: "../../envs/qc.yaml"
    log: f"{OUT}/logs/fastqc/{{sample}}_raw.log"
    shell:
        """
        fastqc --threads {threads} \
               --outdir {params.outdir} \
               {input.r1} {input.r2} \
               > {log} 2>&1
        """


rule fastp_trim:
    """Trim adapters, low-quality bases, and short reads with fastp."""
    input:
        r1 = lambda wc: f"{config['reads_dir']}/{wc.sample}_R1.fastq.gz",
        r2 = lambda wc: f"{config['reads_dir']}/{wc.sample}_R2.fastq.gz",
    output:
        r1      = f"{OUT}/qc/fastp/{{sample}}/{{sample}}_R1_trimmed.fastq.gz",
        r2      = f"{OUT}/qc/fastp/{{sample}}/{{sample}}_R2_trimmed.fastq.gz",
        html    = f"{OUT}/qc/fastp/{{sample}}/{{sample}}_fastp.html",
        json    = f"{OUT}/qc/fastp/{{sample}}/{{sample}}_fastp.json",
    params:
        min_len = config["fastp"]["min_length"],
        qual    = config["fastp"]["qualified_quality_phred"],
        cut_f   = "--cut_front"  if config["fastp"]["cut_front"]  else "",
        cut_t   = "--cut_tail"   if config["fastp"]["cut_tail"]   else "",
        extra   = config["fastp"]["extra_args"],
    threads: 4
    conda: "../../envs/qc.yaml"
    log: f"{OUT}/logs/fastp/{{sample}}.log"
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.qual} \
            {params.cut_f} {params.cut_t} \
            --thread {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule fastqc_trimmed:
    """Run FastQC on trimmed reads (output collected by MultiQC)."""
    input:
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2,
    output:
        html_r1 = f"{OUT}/qc/fastqc/{{sample}}_R1_trimmed_fastqc.html",
        html_r2 = f"{OUT}/qc/fastqc/{{sample}}_R2_trimmed_fastqc.html",
        zip_r1  = f"{OUT}/qc/fastqc/{{sample}}_R1_trimmed_fastqc.zip",
        zip_r2  = f"{OUT}/qc/fastqc/{{sample}}_R2_trimmed_fastqc.zip",
    params:
        outdir = f"{OUT}/qc/fastqc",
    threads: 2
    conda: "../../envs/qc.yaml"
    log: f"{OUT}/logs/fastqc/{{sample}}_trimmed.log"
    shell:
        """
        fastqc --threads {threads} \
               --outdir {params.outdir} \
               {input.r1} {input.r2} \
               > {log} 2>&1
        """
