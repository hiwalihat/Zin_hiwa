# =============================================================================
# Rule module: De Novo Assembly
# Tool: MEGAHIT
# =============================================================================

rule megahit_assembly:
    """Assemble clean reads into contigs with MEGAHIT."""
    input:
        r1 = rules.bowtie2_host_removal.output.r1_clean,
        r2 = rules.bowtie2_host_removal.output.r2_clean,
    output:
        contigs = f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    params:
        outdir       = f"{OUT}/assembly/{{sample}}",
        min_contig   = config["megahit"]["min_contig_len"],
        k_list       = config["megahit"]["k_list"],
        extra        = config["megahit"]["extra_args"],
    threads: config["megahit"]["threads"]
    conda: "../../envs/assembly.yaml"
    log: f"{OUT}/logs/assembly/{{sample}}.log"
    shell:
        """
        # MEGAHIT refuses to write to an existing directory
        rm -rf {params.outdir}

        megahit \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            --min-contig-len {params.min_contig} \
            --k-list {params.k_list} \
            --num-cpu-threads {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule index_contigs:
    """Build a Bowtie2 index of the assembled contigs (needed for binning)."""
    input:
        contigs = rules.megahit_assembly.output.contigs,
    output:
        index_flag = f"{OUT}/assembly/{{sample}}/bowtie2_index/contigs.1.bt2",
    params:
        prefix = f"{OUT}/assembly/{{sample}}/bowtie2_index/contigs",
    threads: 4
    conda: "../../envs/host_removal.yaml"
    log: f"{OUT}/logs/assembly/{{sample}}_bowtie2_index.log"
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {params.prefix} \
            > {log} 2>&1
        """


rule map_reads_to_contigs:
    """Map clean reads back to contigs for coverage estimation (binning)."""
    input:
        r1    = rules.bowtie2_host_removal.output.r1_clean,
        r2    = rules.bowtie2_host_removal.output.r2_clean,
        index = rules.index_contigs.output.index_flag,
    output:
        bam = f"{OUT}/assembly/{{sample}}/reads_vs_contigs.sorted.bam",
        bai = f"{OUT}/assembly/{{sample}}/reads_vs_contigs.sorted.bam.bai",
    params:
        prefix = f"{OUT}/assembly/{{sample}}/bowtie2_index/contigs",
    threads: 8
    conda: "../../envs/host_removal.yaml"
    log: f"{OUT}/logs/assembly/{{sample}}_mapping.log"
    shell:
        """
        bowtie2 \
            -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            2> {log} | \
        samtools sort \
            -@ {threads} \
            -o {output.bam}

        samtools index {output.bam}
        """
