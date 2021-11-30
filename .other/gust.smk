import os
import glob
from pathlib import Path
configfile: "config.yml"
ref_genome = config["reference_genome"]
fragsize = str(config["fragment_size"]) + "bp_fragments/"


# locate input fasta files
fastanames = [Path(i).stem for i in (glob.glob("genomes/*.fasta") + glob.glob("genomes/*.fa"))]
fastagznames = [Path(os.path.splitext(i)[0]).stem for i in (glob.glob("genomes/*.fasta.gz") + glob.glob("genomes/*.fa.gz"))]
allfastanames = fastanames + fastagznames
# not yet implemented
#fastqnames = list(set([Path(i).stem for i in (glob.glob("reads/*.fastq") + glob.glob("reads/*.fq"))]))
#fastqgznames = list(set([Path(os.path.splitext(i)[0]).stem for i in (glob.glob("reads/*.fastq.gz") + glob.glob("reads/*.fq.gz"))]))
#allfastqnames = fastqnames + fastqgznames
#alldatanames = allfastanames + allfastqnames

rule all:
    input:
        variants = fragsize + "snp_discovery/snps.raw.bcf",
        fitlered_variants = fragsize + "snp_discovery/snps.filt.4.bcf"

rule fasta2fastq:
    input:  "genomes/{assembly}.fasta"
    output: "genomes/fastq/{assembly}.fq"
    message: "Using seqtk to convert {input} to FASTQ with quality score J"
    threads: 1
    shell: "seqtk seq -C -U -S -F 'J' {input} > {output}"

rule fragment_assemblies:
    input: "genomes/fastq/{assembly}.fq"
    output: fragsize + "genomes_fragmented/{assembly}.frag.fq.gz"
    params: config["fragment_size"]
    message: "Using seqkit to create {params}bp sliding window fragments from {input}"
    threads: 2
    shell: "seqkit sliding -j {threads} -s 1 -W {params} {input} -o {output}"

rule isolate_reference:
    input: "genomes/" + ref_genome
    output: "genomes/reference/" + ref_genome
    message: "Symlinking {input} to {output}"
    threads: 1
    shell: "ln -sr {input} {output}"

rule index_reference_bwa:
    input: "genomes/reference/" + ref_genome
    output: "genomes/reference/" + ref_genome + ".bwt"
    log: "genomes/reference/" + ref_genome + ".idx.log"
    message: "Indexing {input} with bwa-mem"
    shell:
        """
        bwa index {input} > {log} 2>&1
        """

rule index_reference_samtools:
    input: "genomes/reference/" + ref_genome
    output: "genomes/reference/" + ref_genome + ".fai"
    message: "Indexing {input} with samtools"
    shell:
        """
        samtools faidx {input} > {output}
        """

rule map_frags_to_reference:
    input:
        reference = "genomes/reference/" + ref_genome,
	    idx = "genomes/reference/" + ref_genome + ".bwt",
        query = fragsize + "genomes_fragmented/{assembly}.frag.fq.gz"
    output: temp(fragsize + "alignments/{assembly}.raw.bam")
    log: fragsize + "alignments/logs/{assembly}.log"
    message: "Using bwa to map {input.query} onto {input.reference}"
    threads: 30
    priority: 1
    params: 
        bwa = config["bwa_parameters"],
        header = "{assembly}"
    shell:
        """
        if [ "{threads}" = "2" ]; then
            MAPT=1
            SAMT=1
        else
            MAPT=$(awk "BEGIN {{print {threads}-int({threads}/3)}}")
            SAMT=$(awk "BEGIN {{print int({threads}/3)}}")
        fi
        bwa mem -t $MAPT {params.bwa} -R "@RG\\tID:{params.header}\\tSM:{params.header}\\tPL:Illumina" -a {input.reference} {input.query} 2> {log} | samtools view -@$SAMT -F 0x04 -bh -o {output} -
        """

rule sort_index_alignments:
    input: fragsize + "alignments/{assembly}.raw.bam"
    output:
        bam = fragsize + "alignments/{assembly}.bam",
        idx = fragsize + "alignments/{assembly}.bam.bai"
    message: "Using samtools to sort and index {input}"
    threads: 5
    params: fragsize + "alignments/{assembly}"
    shell:
        """
        samtools sort -@{threads} -T {params} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule alignment_list:
    input: expand(fragsize + "alignments/{assembly}.bam", assembly = allfastanames)
    output: fragsize + "alignments/alignments.list",
    message: "Generating alignment list for variant calling"
    threads: 1
    params: fragsize
    shell:
        """
        ls {params}alignments/*.bam > {output}
        """


#rule merge_alignments:
#    input: expand(fragsize + "alignments/{assembly}.bam", assembly = allfastanames)
#    output:
#        bamlist= fragsize + "alignments/.bamlist",
#        bam = fragsize + "alignments/alignments.bam"
#    message: "Merging all of the alignments of {output.bamlist} into {output.bam}"
#    threads: 30
#    params: fragsize
#    shell:
#        """
#        ls {params}alignments/*.bam > {output.bamlist}
#        samtools merge -@{threads} -b {output.bamlist} -f {output.bam}
#        """
#
#rule index_alignments:
#    input: fragsize + "alignments/alignments.bam"
#    output: fragsize + "alignments/alignments.bam.bai"
#    message: "Indexing {input} with samtools"
#    threads: 1
#    shell: "samtools index {input}"
#
rule create_popmap:
    input: fragsize + "alignments/alignments.list"
    output: fragsize + "populations.map"
    message: "Creating population mapping file {output} based on genome names"
    threads: 1
    shell:
        """
        sed 's/.bam//g' {input} | sed "s/.*\///" | paste -d' ' - - > {output}
        """

rule split_regions:
    input:
        fai =  "genomes/reference/" + ref_genome + ".fai"
    output: fragsize + "snp_discovery/reference.5kb.regions"
    message: "Splitting {input.fai} into 5kb regions for variant calling parallelization"
    threads: 1
    shell:
        """
        fasta_generate_regions.py {input.fai} 5000 > {output}
        """

rule call_variants:
    input:
        alignments = fragsize + "alignments/alignments.list",
        genome = "genomes/reference/" + ref_genome,
        regions = fragsize + "snp_discovery/reference.5kb.regions",
        populations = fragsize + "populations.map"
    output: fragsize + "snp_discovery/snps.raw.bcf"
    message: "Calling variants with freebayes"
    conda: "variantcalling.yml"
    threads: 30
    params: config["freebayes_parameters"]
    shell:
        """
        cat {input.regions} | parallel -k -j {threads} freebayes -f {input.genome} --bam-list {input.alignments} \
            -C 2 --min-coverage 5 --ploidy 1 --standard-filters --populations {input.populations} --region {{}} \
            | vcffirstheader \
            | vcfstreamsort -w 1000 \
            | vcfuniq \
            | bcftools view -Ob - > {output}
        """

rule variant_filter_qual:
    input: fragsize + "snp_discovery/snps.raw.bcf"
    output: fragsize + "snp_discovery/snps.filt.1.bcf"
    message: "Filtering {input} based on genotype quality, mapping quality, and depth"
    shell:
        """
        bcftools view -i'QUAL>=30 && SRR<3 && MQM>40.0 && MQMR>-5.0 && MIN(INFO/DP)>10' {input} > {output}
        """

rule variant_filter_sitedepth:
    input: fragsize + "snp_discovery/snps.filt.1.bcf"
    output: fragsize + "snp_discovery/snps.filt.2.bcf"
    log: fragsize + "snp_discovery/sitedepth.txt"
    message: "Filtering {input} by maximum depth, monomorphism, maximum depth. Also removing sites with any missing data"
    shell:
        """
        bcftools query {input} -f '%DP\n' > {log}
        MAXDP=$(awk '{{ sum += $1; n++ }} END {{ if (n > 0) print sum / n; }}' {log})
        bcftools view -e "F_MISSING > 0.0 || AC==0 || AC=AN || INFO/DP>$MAXDP || INFO/DP < 10" {input} > {output}
        """

rule variant_filter_splitmnp:
    input: fragsize + "snp_discovery/snps.filt.2.bcf"
    output: fragsize + "snp_discovery/snps.filt.3.bcf"
    message: "Splitting SNPs from indels in {input}"
    shell:
        """
        bcftools norm -m -any {input} > {output}
        """

rule variant_genotyping_error:
    input: fragsize + "snp_discovery/snps.filt.3.bcf"
    output: fragsize + "snp_discovery/snps.filt.4.bcf"
    message: "Removing sites where the reference genome self-alignment does not have the reference allele (genotyping error)"
    params: ""
    shell:
        """
        IDX=$(bcftools query -l {input} | awk "/$refname/ {print NR - 1}")
        bcftools filter -s GENOERROR -m + -i "'GT[$IDX]="R"'" {input} > {output}
        """

# needs a way to remove sites with genotyping error where reference genome is not homozygous ref allele for that site
#1 get sample name order in the bcf file
#2 get 0-based index of the reference sample


#3 apply a filter something like this
# where * is the index of the reference sample, which you want to be homozygous for the reference allele

rule variant_filter_LDthinning:
    input: fragsize + "snp_discovery/snps.filt.3.bcf"
    output: fragsize + "snp_discovery/snps.filt.4.bcf"
    params: config["window_size"]
    message: "Thinning SNPs in {input} to retain 1 in every {params}bp"
    shell:
        """
        bcftools +prune -w {params}bp -n 1 -N maxAF {input} > {output}
        """
