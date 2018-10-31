# Parsing the config files
configfile: "config.json"

# Parse the required files
BAMDIR = config["resources"]["bam"]
REF = config["resources"]["reference"]

# Parse path for tools
GATK = config["tools"]["GATK"]
BEAGLE = config["tools"]["BEAGLE"]
BCFTOOLS = config["tools"]["BCFTOOLS"]
LOAD_JAVA = config["tools"]["JAVA"]


# Parse the output directory
OUTDIR = config["base_out"]["out_dir"]

# Parsing wildcards
chr_list = list(range(1, 30))
ch_id = ["Chr" + str(i) for i in chr_list]
SAMPLES, = glob_wildcards(BAMDIR + "/{sample}_recalibrated.bam")

rule all:
    input:
        expand(OUTDIR + "/vcf_raw_beagle/{chromosome}_gatk4_raw_beagle.vcf.gz", chromosome=ch_id),
        expand(OUTDIR + "/vcf_filter_beagle/{chromosome}_gatk4_filter_beagle.vcf.gz", chromosome=ch_id)


# Scatter by sample and chromosome

rule Haplotype_caller:
    input:
        BAMDIR + "/{sample}_recalibrated.bam"
    output:
        OUTDIR + "/gvcf/{chromosome}/{sample}_{chromosome}.g.vcf.gz"
    params:
        chr = "Chr{chromosome}",
        ref = REF
    shell:
        GATK4 +
        " HaplotypeCaller " +
        " -I {input} " +
        " -R {params.ref} " +
        " -L {params.chr} " +
        " -O {output} " +
        " --ERC GVCF "

# Gather all samples given chromosome

rule GenomicsDB_import:
    input:
        expand(OUTDIR + "/gvcf/{{chromosome}}/{sample}_{{chromosome}}.g.vcf.gz", sample=SAMPLES)
    output:
        samp_map = OUTDIR + "/db_imp/{chromosome}.map",
        db_dir = temp(directory(OUTDIR + "/db_imp/db_{chromosome}"))
    params:
        chr = "{chromosome}",
        ref = REF
    shell:
        " paste -- " +
        " <(echo {input}|xargs -n 1 basename -s \"_{params.chr}.g.vcf.gz\") " +
        " <(echo {input}|tr ' ' '\\n' ) > {output.samp_map} \n " +
        GATK4 +
        " GenomicsDBImport " +
        " --sample-name-map {output.samp_map} " +
        " --genomicsdb-workspace-path {output.db_dir} " +
        " -L {params.chr} " +
        " -R {params.ref} "

rule Joint_Genotyping:
    input:
        OUTDIR + "/db_imp/db_{chromosome}"
    output:
        OUTDIR + "/vcf/{chromosome}_gatk4.vcf.gz"
    params:
        chr = "{chromosome}",
        ref = REF
    shell:
        GATK4 +
        " GenotypeGVCFs " +
        " -R {params.ref} " +
        " -L {params.chr} " +
        " -O {output} " +
        " -V gendb://{input} "

rule raw_beagle_imputation:
    input:
        OUTDIR + "/vcf/{chromosome}_gatk4.vcf.gz"
    output:
        OUTDIR + "/vcf_raw_beagle/{chromosome}_gatk4_raw_beagle.vcf.gz"
    shell:
        " java -jar " + BEAGLE +
        " gl={input} " +
        " out=" + OUTDIR + "/vcf_beagle/{wildcards.chromosome}_gatk4_beagle "


rule select_SNP:
    input:
        OUTDIR + "/vcf/{chromosome}_gatk4.vcf.gz"
    output:
        temp(OUTDIR + "/snp/{chromosome}_gatk4_snp.vcf.gz")
    params:
        ref = REF
    shell:
        GATK4 +
        " SelectVariants  " +
        " -R {params.ref} " +
        " -V {input}  " +
        " --select-type-to-include SNP  " +
        " --output {output} "

rule filter_SNP:
    input:
        rules.select_SNP.output
    output:
        temp(OUTDIR + "/snp/{chromosome}_gatk4_snp_filtered.vcf.gz")
    params:
        ref = REF
    shell:
        GATK4 +
        "  VariantFiltration " +
        "  -R {params.ref} " +
        "  -V {input} " +
        "  --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" " +
        "  --filter-name \"snp_filter\" " +
        "  --output {output} "

rule select_indel:
    input:
        OUTDIR + "/vcf/{chromosome}_gatk4.vcf.gz"
    output:
        temp(OUTDIR + "/indel/{chromosome}_gatk4_indel.vcf.gz")
    params:
        ref = REF
    shell:
        GATK4 +
        " SelectVariants  " +
        " -R {params.ref} " +
        " -V {input}  " +
        " --select-type-to-include INDEL " +
        " --output {output} "


rule filter_indel:
    input:
        rules.select_indel.output
    output:
        temp(OUTDIR + "/indel/{chromosome}_gatk4_indel_filtered.vcf.gz")
    params:
        ref = REF
    shell:
        GATK4 +
        "  VariantFiltration " +
        "  -R {params.ref} " +
        "  -V {input} " +
        "  --filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\" " +
        "  --filter-name \"indel_filter\" " +
        "  --output {output} "


rule merge_variants:
    input:
        sic = rules.filter_SNP.output,
        idc = rules.filter_indel.output
    output:
        OUTDIR + "/vcf_filtered/{chromosome}_gatk4_mergefil.vcf.gz"
    shell:
        GATK4 +
        " MergeVcfs " +
        " --INPUT {input.sic} " +
        " --INPUT {input.idc} " +
        " --OUTPUT {output} "

rule remove_filtered:
    input:
        rules.merge_variants.output
    output:
        OUTDIR + "/vcf_filtered/{chromosome}_gatk4_pruned.vcf.gz"
    params:
        ref = REF
    shell:
        BCFTOOLS +
        " view " +
        " -f PASS " +
        " -o {output} -O z " +
        " {input} "

rule filtered_beagle_imputation:
    input:
        rules.remove_filtered.output
    output:
        OUTDIR + "/vcf_filter_beagle/{chromosome}_gatk4_filter_beagle.vcf.gz"
    shell:
        " java -jar " + BEAGLE +
        " gl={input} " +
        " out=" + OUTDIR + "/vcf_filter_beagle/{wildcards.chromosome}_gatk4_filter_beagle"
