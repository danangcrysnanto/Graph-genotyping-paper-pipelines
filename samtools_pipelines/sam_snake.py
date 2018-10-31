# Load the required python library
configfile: "config.json"
import vcf
import statistics

# Parse the required files
BAMDIR = config["resources"]["bam"]
REF = config["resources"]["reference"]

# Parse tools
SAMTOOLS = config["tools"]["SAMTOOLS"]
BCFTOOLS = config["tools"]["BCFTOOLS"]
GATK = config["tools"]["GATK"]
BEAGLE = config["tools"]["BEAGLE"]
BGZIP = config["tools"]["BGZIP"]
TABIX = config["tools"]["TABIX"]
LOAD_JAVA = config["tools"]["JAVA"]


# Parsing output directory
OUTDIR = config["base_out"]["out_dir"]

# Parsing wildcards
ch_list = list(range(1,30))
SAMPLES, = glob_wildcards(BAMDIR + "/{sample}_recalibrated.bam")


rule all:
    input:
        expand(OUTDIR + "/vcf_filtered_beagle/{chromosome}_samtools_filtered_beagle.vcf.gz", chromosome=ch_list),
        expand(OUTDIR + "/vcf_raw_beagle/{chromosome}_samtools_raw_beagle.vcf.gz", chromosome=ch_list)


rule calling_and_genotype:
    input:
        expand(BAMDIR + "{sample}_recalibrated.bam", sample=SAMPLES)
    output:
        OUTDIR + "/vcf_raw/{chromosome}_samtools.vcf.gz"
    shell:
        SAMTOOLS +  " mpileup -uf " + REF +
        " --region Chr{wildcards.chromosome} " +
        " -E " +
        " -t DP -t SP -t ADF -t ADR -t AD {input}| " +
        BCFTOOLS + " call -mv -o {output} -O z"


rule overlap_filter:
    input:
        OUTDIR + "/vcf_raw/{chromosome}_samtools.vcf.gz"
    output:
        temp(OUTDIR + "/vcf_nooverlap/{chromosome}_overlap.txt"),
        temp(OUTDIR + "/vcf_nooverlap/{chromosome}_samtools_nooverlap.vcf.gz")
    shell:
        "awk \' ! /^#/{{print $2}} \' <(zcat {input})| " +
        "sort -n| " +
        "uniq -d > {output[0]} \n"
        "awk \' NR==FNR{{arr[$0];next}}/^#/ {{print $0}} ! ($2 in arr){{print $0}} \' {output[0]} <(zcat {input})| " + 
        "bgzip > {output[1]} && tabix {output[1]}"


def get_depth_threshold(input_file):
    vcf_reader = vcf.Reader(open(input_file, 'rb'))
    dp_list = []
    for record in vcf_reader:
            dp_list.append(record.INFO["DP"])
    return statistics.median(dp_list) + (3 * statistics.mean(dp_list))


rule quality_filter:
    input:
        OUTDIR + "/vcf_nooverlap/{chromosome}_samtools_nooverlap.vcf.gz"
    output:
        OUTDIR + "/vcf_filtered/{chromosome}_samtools_filtered.vcf.gz"
    run:
        depth_limit = get_depth_threshold(input[0])
        shell(GATK + " VariantFiltration " +
              " -R " + REF +
              " -V {input} " +
              " --filterExpression \"QUAL < 20 || MQ < 30 || DP < 10 || DP > " + str(depth_limit) + " \"" +
              " --filterName \"1000_bull_filter \" " +
              " -o {output} ")


rule beagle_imputation_raw:
    input:
        OUTDIR + "/vcf_raw/{chromosome}_samtools.vcf.gz"
    output:
        OUTDIR + "/vcf_raw_beagle/{chromosome}_samtools_raw_beagle.vcf.gz"
    shell:
        " java -jar " + BEAGLE +
        " gl={input} " +
        " out=" + OUTDIR + "/vcf_beagle/{wildcards.chromosome}_samtools_raw_beagle "


rule remove_filtered:
    input:
        OUTDIR + "/vcf_filtered/{chromosome}_samtools_filtered.vcf.gz"
    output:
        OUTDIR + "/vcf_filtered_pruned/{chromosome}_samtools_pruned.vcf.gz"
    shell:
        BCFTOOLS +
        " view " +
        " -f PASS " +
        " -o {output} -O z " +
        " {input} "


rule beagle_imputation_filtered:
    input:
        OUTDIR + "/vcf_filtered_pruned/{chromosome}_samtools_pruned.vcf.gz"
    output:
        OUTDIR + "/vcf_filtered_beagle/{chromosome}_samtools_filtered_beagle.vcf.gz"
    shell:
        " java -jar " + BEAGLE +
        " gl={input} " +
        " out=" + OUTDIR + "/vcf_filtered_beagle/{wildcards.chromosome}_samtools_filtered_beagle "