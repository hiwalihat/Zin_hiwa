# =============================================================================
# Rule module: Functional Profiling
# Tool: HUMAnN3
# =============================================================================

rule humann3_profile:
    """
    Run HUMAnN3 to produce pathway abundances, gene family abundances,
    and pathway coverage tables.
    """
    input:
        # HUMAnN3 can work on a single (concatenated) FASTQ; concatenate R1+R2
        r1 = rules.bowtie2_host_removal.output.r1_clean,
        r2 = rules.bowtie2_host_removal.output.r2_clean,
    output:
        genefamilies    = f"{OUT}/functional/{{sample}}/{{sample}}_genefamilies.tsv",
        pathabundance   = f"{OUT}/functional/{{sample}}/{{sample}}_pathabundance.tsv",
        pathcoverage    = f"{OUT}/functional/{{sample}}/{{sample}}_pathcoverage.tsv",
    params:
        outdir          = f"{OUT}/functional/{{sample}}",
        nuc_db          = config["databases"]["humann3_nucleotide_db"],
        prot_db         = config["databases"]["humann3_protein_db"],
        metaphlan_opts  = config["humann3"]["metaphlan_options"],
        extra           = config["humann3"]["extra_args"],
        cat_reads       = f"{OUT}/functional/{{sample}}/{{sample}}_cat.fastq.gz",
    threads: config["humann3"]["threads"]
    conda: "../../envs/functional.yaml"
    log: f"{OUT}/logs/functional/{{sample}}.log"
    shell:
        """
        # Concatenate paired reads for HUMAnN3 input
        cat {input.r1} {input.r2} > {params.cat_reads}

        humann \
            --input {params.cat_reads} \
            --output {params.outdir} \
            --output-basename {wildcards.sample} \
            --nucleotide-database {params.nuc_db} \
            --protein-database {params.prot_db} \
            --metaphlan-options "{params.metaphlan_opts}" \
            --threads {threads} \
            {params.extra} \
            > {log} 2>&1

        # Remove large intermediate concatenated file
        rm -f {params.cat_reads}
        """


rule humann3_renorm:
    """Renormalise gene families and pathway abundances to relative abundance (copies-per-million)."""
    input:
        genefamilies  = rules.humann3_profile.output.genefamilies,
        pathabundance = rules.humann3_profile.output.pathabundance,
    output:
        genefamilies_cpm  = f"{OUT}/functional/{{sample}}/{{sample}}_genefamilies_cpm.tsv",
        pathabundance_cpm = f"{OUT}/functional/{{sample}}/{{sample}}_pathabundance_cpm.tsv",
    conda: "../../envs/functional.yaml"
    log: f"{OUT}/logs/functional/{{sample}}_renorm.log"
    shell:
        """
        humann_renorm_table \
            --input {input.genefamilies} \
            --output {output.genefamilies_cpm} \
            --units cpm \
            >> {log} 2>&1

        humann_renorm_table \
            --input {input.pathabundance} \
            --output {output.pathabundance_cpm} \
            --units cpm \
            >> {log} 2>&1
        """
