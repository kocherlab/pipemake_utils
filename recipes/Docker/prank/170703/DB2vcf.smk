chrom_list = ['MGC2_chr_1', 'MGC2_chr_2', 'MGC2_chr_3', 'MGC2_chr_6', 'MGC2_chr_8', 'MGC2_chr_9']

rule all:
    input:
        expand("results/MGC2_v3/vcfs/raw_{chrom}.vcf.gz", chrom=chrom_list),
        expand("results/MGC2_v3/vcfs/raw_{chrom}.vcf.gz.tbi", chrom=chrom_list)

rule DB2vcf:
    """
    This rule uses the genomic databases from the previous step (gvcf2DB) to create VCF files, one per list file. Thus, lists
    are still scattered.
    """
    input:
        db = "results/MGC2_v3/genomics_db_import/{chrom}_DB.tar",
        ref = "results/MGC2_v3/data/genome/MGC2_v3.fna",
        intervals = "results/MGC2_v3/genomics_db_import/intervals/{chrom}.list"
    output:
        vcf = temp("results/MGC2_v3/vcfs/raw_{chrom}.vcf.gz"),
        vcfidx = temp("results/MGC2_v3/vcfs/raw_{chrom}.vcf.gz.tbi")
    params:
        het = config['het_prior'],
        db = lambda wc, input: input.db[:-4]
    resources:
        mem_mb = 128000
        reduced = 125000
    log:
        "logs/MGC2_v3/gatk_genotype_gvcfs/{chrom}.txt"
    benchmark:
        "benchmarks/MGC2_v3/gatk_genotype_gvcfs/{chrom}.txt"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        """
        tar -xf {input.db}
        gatk GenotypeGVCFs \
            --java-options '-Xmx{resources.reduced}m -Xms{resources.reduced}m' \
            -R {input.ref} \
            --heterozygosity {params.het} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --all-sites \
            -L {input.intervals} \
            --tmp-dir {resources.tmpdir} &> {log}
        """