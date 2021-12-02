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
        fitlered_variants = fragsize + "snp_discovery/snps.filt.5.bcf",
        msa_input = fragsize + "msa/filtered.variants.fasta",
        msa = fragsize + "msa/variants.resampled.efa"

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
        bwa mem -t $MAPT {params.bwa} -R "@RG\\tID:{params.header}\\tSM:{params.header}\\tPL:Illumina" {input.reference} {input.query} 2> {log} | samtools view -@$SAMT -F 0x04 -bh -o {output} -
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
        echo {input} | tr " " "\n" > {output}
        """

rule create_popmap:
    input: fragsize + "alignments/alignments.list"
    output: fragsize + "populations.map"
    message: "Creating population mapping file {output} based on genome names"
    threads: 1
    shell:
        """
        paste <(sed 's/.bam//g' {input} | sed "s/.*\///") <(sed 's/.bam//g' {input} | sed "s/.*\///") > {output}
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
    log: fragsize + "snp_discovery/variant.stats"
    message: "Calling variants with freebayes"
    threads: 30
    params: config["freebayes_parameters"]
    shell:
        """
        cat {input.regions} | parallel -k -j {threads} freebayes -f {input.genome} --bam-list {input.alignments} \
            {params} --populations {input.populations} --region {{}} \
            | vcffirstheader \
            | vcfstreamsort -w 1000 \
            | vcfuniq \
            | bcftools view -Ob - > {output}
        echo "file,sites,filter" > {log}
        echo -n "$(basename {output} .bcf),$(bcftools query {output} -f '%DP\n' | wc -l)" >> {log}
        echo ",raw variant sites" >> {log}
        """

rule variant_filter_qual:
    input: fragsize + "snp_discovery/snps.raw.bcf"
    output: fragsize + "snp_discovery/snps.filt.1.bcf"
    message: "Filtering {input} based on genotype quality, mapping quality, and depth"
    params: fragsize
    shell:
        """
        bcftools view -i'QUAL>=30 && SRR<3 && MQM>40.0 && MQMR>-5.0 && MIN(INFO/DP)>10' {input} > {output}
        echo -n "$(basename {output} .bcf),$(bcftools query {output} -f '%DP\n' | wc -l)" >> {params}snp_discovery/variant.stats
        echo ",include QUAL>=30 && SRR<3 && MQM>40.0 && MQMR>-5.0 && MIN(INFO/DP)>10" >> {params}snp_discovery/variant.stats
        """

rule variant_filter_sitedepth:
    input: fragsize + "snp_discovery/snps.filt.1.bcf"
    output: fragsize + "snp_discovery/snps.filt.2.bcf"
    log: fragsize + "snp_discovery/sitedepth.txt"
    message: "Filtering {input} by maximum depth, monomorphism, maximum depth. Also removing sites with any missing data"
    params: fragsize
    shell:
        """
        bcftools query {input} -f '%DP\n' > {log}
        MAXDP=$(awk '{{ sum += $1; n++ }} END {{ if (n > 0) print sum / n; }}' {log})
        bcftools view -e "F_MISSING > 0.0 || AC==0 || AC=AN || INFO/DP>$MAXDP || INFO/DP < 10" {input} > {output}
        echo -n "$(basename {output} .bcf),$(bcftools query {output} -f '%DP\n' | wc -l)" >> {params}snp_discovery/variant.stats
        echo ",exclude F_MISSING > 0.0 || AC==0 || AC=AN || INFO/DP>$MAXDP || INFO/DP < 10" >> {params}snp_discovery/variant.stats
        """

rule variant_filter_splitmnp:
    input: fragsize + "snp_discovery/snps.filt.2.bcf"
    output: fragsize + "snp_discovery/snps.filt.3.bcf"
    message: "Splitting SNPs from indels in {input} and removing indels"
    params: fragsize
    shell:
        """
        bcftools norm -m -any -a {input} | bcftools view --types snps > {output}
        echo -n "$(basename {output} .bcf),$(bcftools query {output} -f '%DP\n' | wc -l)" >> {params}snp_discovery/variant.stats
        echo ",norm -m -any -a and --types snps" >> {params}snp_discovery/variant.stats
        """

# needs a way to remove sites with genotyping error where reference genome is not homozygous ref allele for that site
#1 get sample name order in the bcf file
#2 get 0-based index of the reference sample
#3 apply a filter something like this
# where * is the index of the reference sample, which you want to be homozygous for the reference allele
rule variant_genotyping_error:
    input: fragsize + "snp_discovery/snps.filt.3.bcf"
    output: fragsize + "snp_discovery/snps.filt.4.bcf"
    message: "Removing sites where the reference genome self-alignment does not have the reference allele (genotyping error)"
    params: 
        refgeno = ref_genome,
        dir = fragsize
    shell:
        """
        refname=$(basename {params.refgeno} .fasta)
        IDX=$(bcftools query -l {input} | awk "/$refname/ {{print NR - 1}}")
        bcftools filter -s GENOERROR -m + -i "'GT[$IDX]="R"'" {input} > {output}
        echo -n "$(basename {output} .bcf),$(bcftools query {output} -f '%DP\n' | wc -l)" >> {params.dir}snp_discovery/variant.stats
        echo ",include GT[$IDX]=\"R\"" >> {params.dir}snp_discovery/variant.stats
        """


rule variant_filter_LDthinning:
    input: fragsize + "snp_discovery/snps.filt.4.bcf"
    output: fragsize + "snp_discovery/snps.filt.5.bcf"
    params: 
        window = config["window_size"],
        dir = fragsize
    message: "Thinning SNPs in {input} to retain 1 in every {params.window}bp"
    shell:
        """
        bcftools +prune -w {params.window}bp -n 1 -N maxAF {input} > {output}
        echo -n "$(basename {output} .bcf),$(bcftools query {output} -f '%DP\n' | wc -l)" >> {params.dir}/snp_discovery/variant.stats
        echo ",+prune -w {params.window}bp -n 1 -N maxAF" >> {params.dir}/snp_discovery/variant.stats
        # nicer fixed-width table
        column -t -s"," {params.dir}/snp_discovery/variant.stats > {params.dir}snp_discovery/.variant.stats \
        && rm {params.dir}/snp_discovery/variant.stats \
        && mv {params.dir}/snp_discovery/.variant.stats {params.dir}snp_discovery/variant.stats
        """

rule vcf2fasta:
    input: fragsize + "snp_discovery/snps.filt.5.bcf"
    output: fragsize + "msa/filtered.variants.fasta"
    message: "Converting the filtered variants into FASTA format for multiple sequence alignment (MSA)"
    params:
        outgroup = config["outgroup"],
        outdir = fragsize + "msa"
    threads: 1
    shell: 
        """
        outgroup=$(basename {params.outgroup} .fasta)
        tools/vcf2phylip.py -i {input} -p -f -m 2 --output-folder {params.outdir} --output-prefix filtered.variants -o $outgroup
        mv {params.outdir}/filtered.variants.min2.fasta {params.outdir}/filtered.variants.fasta
        """

rule muscle_msa:
    input: fragsize + "msa/filtered.variants.fasta"
    output: fragsize + "msa/variants.diversified.efa"
    log: fragsize + "msa/variants.diversified.log"
    message: "Using MUSCLE to perform diversified multiple sequence alignment (MSA). This will likely take several hours."
    params: config["musclev5_parameters"]
    threads: 30
    shell:
        """
        tools/muscle_v5.0.1428_linux -align {input} -diversified -output {output} -nt -threads {threads} {params}
        tools/muscle_v5.0.1428_linux -disperse {output} -log {log}
        """

rule extract_best_msa:
    input: fragsize + "msa/variants.diversified.efa"
    output: fragsize + "msa/best.msa.afa"
    message: "Extracting best MSA tree"
    threads: 1
    shell: "tools/muscle_v5.0.1428_linux -maxcc {input} -output {output}"

rule resample_msa:
    input: fragsize + "msa/variants.diversified.efa"
    output: fragsize + "msa/variants.resampled.efa" #afa = fragsize + "msa/variants.resampled.afa"
    message: "Resampling diversified MSA ensemble"
    shell: 
        """
        tools/muscle_v5.0.1428_linux -resample {input} -output {output}
        tools/muscle_v5.0.1428_linux -efa_exploide {output}
        """

#rule build_tree_raxml:
#    input: fragsize + "msa/variants.msa.phy"
#    output: fragsize + "phylo/phylogenies.phy"
#    message:
#    threads:
#    params:
#    shell:
#        """
#        raxml-ng -f a -p 12345 -x 12345 -N 20 -s {input} -m GTRCAT -n boot1
#        raxml-ng -f b -t ref -z RAxML_bootstrap.boot1 -m GTRCAT -n consensus
#        """
