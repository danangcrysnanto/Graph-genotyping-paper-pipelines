# Parse config files
configfile: "config.json"

# Parse the required files and parameter

# Parse the required files
GENOME = config["resources"]["reference"]
REF_INDEX = config["resources"]["index_ref"]
BAMDIR = config["resources"]["bam"]


# Parse the tools
GRAPHTYPER = config["tools"]["Graphtyper"]
SAMTOOLS = config["tools"]["Samtools"]
VT = config["tools"]["VT"]
VCFFILTER = config["tools"]["vcffilter"]
BEAGLE = config["tools"]["Beagle"]
BGZIP = config["tools"]["BGZIP"]
TABIX = config["tools"]["tabix"]

# Parse the chunk size and padding
slc = config["run"]["slice"]
padding = config["run"]["padding"]
CALLOPTS = config["run"]["opts"]

chr_list = list(range(1, 30))
samples, = glob_wildcards(BAMDIR + "/{sample}_recalibrated.bam")

rule all:
    input:
        expand("vcf_raw_beagle/vcf_raw_beagle_{chromosome}.vcf.gz", chromosome=chr_list),
        expand("vcf_filtered_beagle/vcf_filtered_beagle_{chromosome}.vcf.gz", chromosome=chr_list)


### Get the padded region for given job id 
def getstartpos(wildcards):
    js = int(wildcards.job_id)
    if js == 1:
        start = 1
    else:
        start = ((js-1)*slc) - padding
    return start


def getstoppos(wildcards):
    js = int(wildcards.job_id)
    with open(REF_INDEX,"r") as bg:
        chr_length=int(bg.readlines()[int(wildcards.chromosome)-1].strip().split('\t')[1])
    if js*slice >= chr_length:
        stop = chr_length
    else:
        stop = (js*slc)+padding
    return stop


# This rule to split the BAM given region 
rule split_bam:
    input:
        "bam/{sample}_recalibrated.bam"
    output:
        "bam_small/{chromosome}/{job_id}/{sample}_{chromosome}_{job_id}.bam"
    params:
        chr = "{chromosome}",
        start = getstartpos,
        stop = getstoppos,
        samtools = SAMTOOLS
    shell:
        """

        {params.samtools} view -bh -o {output} -O z {input} chr{params.chr}:{params.start}-{params.stop}
       
        """
rule get_bam_list:
    input:
        expand("bam_small/{{chromosome}}/{{job_id}}/{sample}_{{chromosome}}_{{job_id}}.bam", sample=samples)
    output:
        "bam_small/{chromosome}/{job_id}/{chromosome}_{job_id}.list"
    params:
        chr = "{chromosome}",
        jobid = "{job_id}"
    shell:
        """
        
        ls {input}  > {output}

        """

rule discovery_linear:
    input:
        rules.get_bam_list.output
    output:
        vcf = "iteration/{chromosome}/{job_id}/D1/{chromosome}_{job_id}_D1.vcf.gz",
        vcf_index = "iteration/{chromosome}/{job_id}/D1/{chromosome}_{job_id}_D1.vcf.gz.tbi"
    params:
        graphtyper = GRAPHTYPER,
        vt = VT,
        tabix = TABIX,
        bgzip = BGZIP,
        dir_jobid = "iteration/{chromosome}/{job_id}/D1",
        chr = "{chromosome}",
        jobid = "{job_id}",
        genome = GENOME,
        start = getstartpos,
        stop = getstoppos,
        opt = CALLOPTS
    shell:
        """

        {params.graphtyper} construct {params.dir_jobid}/graph {params.genome} chr{params.chr}:{params.start}-{params.stop}
       
        {params.graphtyper} index {params.dir_jobid}/graph

        {params.graphtyper} call {params.opt} {params.dir_jobid}/graph "." \
        --output {params.dir_jobid} \
        --sams  {input} 
        
        rm -f {params.dir_jobid}/graph && rm -rf {params.dir_jobid}/graph_gti

        {params.vt} sort -o {params.dir_jobid}/new_region_sorted.vcf.gz {params.dir_jobid}/*_variants.vcf.gz
        {params.vt} uniq -o {params.dir_jobid}/new.vcf  {params.dir_jobid}/new_region_sorted.vcf.gz 
        cat {params.dir_jobid}/new.vcf | {params.bgzip} -c > {output[0]}
        {params.tabix} {output[0]}

        """

rule discovery_graph:
    input:
        rules.get_bam_list.output,
        rules.discovery_linear.output.vcf
    output:
        vcf = "iteration/{chromosome}/{job_id}/D2/{chromosome}_{job_id}_D2.vcf.gz",
        vcf_index = "iteration/{chromosome}/{job_id}/D2/{chromosome}_{job_id}_D2.vcf.gz.tbi"
    params:
        graphtyper = GRAPHTYPER,
        vt = VT,
        tabix = TABIX,
        bgzip= BGZIP,
        dir_jobid = "iteration/{chromosome}/{job_id}/D2",
        chr = "{chromosome}",
        genome = GENOME,
        start = getstartpos,
        stop = getstoppos,
        opt = CALLOPTS
    shell:
        """

        {params.graphtyper} construct {params.dir_jobid}/graph {params.genome} --vcf={input[1]} chr{params.chr}:{params.start}-{params.stop}
        {params.graphtyper} index {params.dir_jobid}/graph

        {params.graphtyper} call {params.opt} {params.dir_jobid}/graph "." \
        --output {params.dir_jobid} \
        --sams  {input[0]} 
        
        
        ls {params.dir_jobid}/*.hap > {params.dir_jobid}/haps2

        {params.graphtyper} haplotypes {params.dir_jobid}/graph \
        --haplotypes {params.dir_jobid}/graph \
        --output {params.dir_jobid}/haps.vcf.gz 

        rm --force {params.dir_jobid}/graph && rm -r --force {params.dir_jobid}/graph_gti
    
        vt cat {params.dir_jobid}/haps.vcf.gz {params.dir_jobid}/*_variants.vcf.gz |
        vt sort -o {params.dir_jobid}/new_region_sorted2.vcf.gz - 
        vt uniq -o {params.dir_jobid}/new2.vcf {params.dir_jobid}/new_region_sorted2.vcf.gz 
        cat {params.dir_jobid}/new2.vcf | {params.bgzip} -c > {output.vcf}
        tabix {output.vcf}
        
        """


rule cleaning_graph:
    input:
        rules.get_bam_list.output,
        rules.discovery_graph.output.vcf
    output:
        vcf = "iteration/{chromosome}/{job_id}/G1/{chromosome}_{job_id}_haps.vcf.gz",
        vcf_index = "iteration/{chromosome}/{job_id}/G1/{chromosome}_{job_id}_haps.vcf.gz.tbi"
    params:
        graphtyper = GRAPHTYPER,
        vt = VT,
        tabix = TABIX,
        dir_jobid = "iteration/{chromosome}/{job_id}/G1",
        chr = "{chromosome}",
        jobid = "{job_id}",
        genome = GENOME,
        start = getstartpos,
        stop = getstoppos,
        opt = CALLOPTS
    shell:
        """ 
        {params.graphtyper} construct {params.dir_jobid}/graph {params.genome} --vcf={input[1]} chr{params.chr}:{params.start}-{params.stop}
        {params.graphtyper} index {params.dir_jobid}/graph
        {params.graphtyper} call {params.opt} {params.dir_jobid}/graph "." \
            --no_new_variants \
            --output {params.dir_jobid} \
            --sams input[0] \

        ls {params.dir_jobid}/*.hap > {params.dir_jobid}/haps3

        {params.graphtyper} haplotypes {params.dir_jobid}/graph \
            --haplotypes {params.dir_jobid}/haps3 \
            --output {output.vcf} \
            --skip_breaking_down_extracted_haplotypes

        tabix {output.vcf}

        rm -f {params.dir_jobid}/graph && rm -rf {params.dir_jobid}/graph_gti/

        """


rule genotyping_variants:
    input:
        rules.get_bam_list.output,
        rules.cleaning_graph.output.vcf
    output:
        vcf = "iteration/{chromosome}/{job_id}/G2/{chromosome}_{job_id}_haps.vcf.gz",
        vcf_index = "iteration/{chromosome}/{job_id}/G2/{chromosome}_{job_id}_haps.vcf.gz.tbi"
    params:
        graphtyper = GRAPHTYPER,
        vt = VT,
        tabix = TABIX,
        bgzip = BGZIP,
        dir_jobid = "iteration/{chromosome}/{job_id}/G2",
        chr = "{chromosome}",
        jobid = "{job_id}",
        genome = GENOME,
        start = getstartpos,
        stop = getstoppos,
        opt = CALLOPTS
    shell:
        """

        {params.graphtyper} construct {params.dir_jobid}/graph {params.genome} --vcf={input[1]} chr{params.chr}:{params.start}-{params.stop}
        {params.graphtyper} index $GRAPH

        {params.graphtyper} call {params.opt} {params.dir_jobid}/graph "." \
             --no_new_variants \
             --output {output.vcf} \
             --sams  {input[0]}

        {params.graphtyper} vcf_break_down \
            {params.dir_jobid}/graph {params.dir_jobid}/*_calls.vcf.gz | 
            bgzip -c > {output.vcf} && tabix {output.vcf}

        """

#Get the padded region for given job id 
def get_nopad_start(wildcards):
    js = int(wildcards.job_id)
    if js == 1:
        start_nopad = 1
    else:
        start_nopad = ((js - 1) * slc) - 1
    return start_nopad

def get_nopad_stop(wildcards):
    js = int(wildcards.job_id)
    with open(REF_INDEX,"r") as bg:
        chr_length=int(bg.readlines()[int(wildcards.chromosome)-1].strip().split('\t')[1])
    if js*slice >= chr_length:
        stop_nopad = chr_length
    else:
        stop_nopad = (js * slc) - 2
    return stop_nopad


rule repadding_variants:
    input:
        rules.genotyping_variants.output.vcf
    output:
       vcf = "iteration/{chromosome}/{job_id}/repad/{chromosome}_{job_id}_repad.vcf.gz",
       vcf_index = "iteration/{chromosome}/{job_id}/repad/{chromosome}_{job_id}_repad.vcf.gz.tbi"
    params:
       chr = "{chromosome}",
       start_nopad = get_nopad_start,
       stop_nopad = get_nopad_stop,
       tabix = TABIX,
       bgzip = BGZIP
    shell:
        """
        {params.tabix} {input} -h chr{params.chr}:{params.start_nopad}-{params.stop_nopad} |
        {params.bgzip} -c > {output.vcf}

        {params.tabix} {output.vcf}

        """

def get_all_input_given_chromosome(wildcards):
    list_file = []
    start = 1
    count = 0
    with open(REF_INDEX,"r") as bg:
        chr_length=int(bg.readlines()[int(wildcards.chromosome)-1].strip().split('\t')[1])
        while start < chr_length:
            stop = start + slice
            job_id = count + 1
            count += 1
            if stop >= chr_length:
                stop = chr_length
            list_file.append("iteration/"+wildcards.chromosome+"/"+str(job_id)+"/repad/"+wildcards.chromosome+"_"+str(job_id)+"_repad.vcf.gz")
            start = stop
    return(list_file)


rule combine_chromosome_level_variants:
    input:
        get_all_input_given_chromosome
    output:
        vcf = "vcf/vcf_united_{chromosome}.vcf.gz",
        vcf_index = "vcf/vcf_united_{chromosome}.vcf.gz.tbi"
    params:
        vt = VT,
        tabix = TABIX,
        bgzip = BGZIP,
        genome = GENOME
    shell:
         """

         {params.vt} cat {input} | 
         {params.vt} sort - | 
         {params.vt} uniq - |
         {params.vt} normalize - |
         {params.bgzip} > {output.vcf} 

         {params.tabix} {output.vcf}

         """ 

rule filtering_variants:
    input:
        rules.combine_chromosome_level_variants.output
    output:
        vcf = "vcf_filtered/vcf_filtered_{chromosome}.vcf.gz",
        vcf_index = "vcf_filtered/vcf_filtered_{chromosome}.vcf.gz.tbi"
    params:
        vcffilter = VCFFILTER,
        bgzip = BGZIP,
        tabix = TABIX
    shell:
        """
        
        {params.vcffilter} \
        -f "ABHet < 0.0 | ABHet > 0.33" -f "ABHom < 0.0 | ABHom > 0.97" -f "MaxAASR > 0.4" -f "MQ > 30"  \
        -F "graphtyper_filter" \
        -t "pass" {input}  |
        {params.bgzip} > {output.vcf}

        {params.tabix} {output.vcf}
        
        """

rule raw_beagle_imputation:
    input:
        rules.combine_chromosome_level_variants.output.vcf
    output:
        vcf = "vcf_raw_beagle/vcf_raw_beagle_{chromosome}.vcf.gz",
        vcf_index = "vcf_raw_beagle/vcf_raw_beagle_{chromosome}.vcf.gz.tbi"
    params:
        chr = "{chromosome}",
        beagle = BEAGLE
    shell:
        """

        java -jar {params.beagle} \
        gl={input} \
        out=vcf_raw_beagle/vcf_raw_beagle_{params.chr}

        """

rule remove_filtered:
    input:
        rules.filtering_variants.output
    output:
        vcf = "vcf_filtered_removed/vcf_filtered_removed_{chromosome}.vcf.gz",
        vcf_index = "vcf_filtered_removed/vcf_filtered_removed_{chromosome}.vcf.gz.tbi"   
    params:
        vt = VT,
        bgzip = BGZIP,
        tabix = TABIX
    shell:
        """
        
        {params.vt} -f "filter.pass" {input} |
        {params.bgzip} > {output.vcf}

        {params.tabix} {output.vcf}

        """

rule filtered_beagle_imputation:
    input:
        rules.remove_filtered.output.vcf
    output:
        vcf = "vcf_filtered_beagle/vcf_filtered_beagle_{chromosome}.vcf.gz",
        vcf_index = "vcf_filtered_beagle/vcf_filtered_beagle_{chromosome}.vcf.gz.tbi" 
    params:
        chr = "{chromosome}",
        beagle = BEAGLE
    shell:
        """

        java -jar {params.beagle} \
        gl={input} \
        out=vcf_filtered_beagle/vcf_filtered_beagle_{params.chr}

        """
