import os
import glob
from pathlib import Path
configfile: "config.yml"
ref_genome = "genomes/" + config["reference_genome"]
fragsize = str(config["fragment_size"]) + "bp_fragments/"

def assemblynames():
    names = glob.glob("genomes/*.fasta")
    gzstrippednames = []
    for i in names:
        if i.endswith(".gz"):
            gzstrippednames.append(os.path.splitext(i)[0])
        else:
            gzstrippednames.append(i)
    return [Path(i).stem for i in gzstrippednames]
assemblies = assemblynames()


rule all:
    input: 
        variants = fragsize + "variants/snps.raw.bcf"

rule fastq_convert:
    input:  "genomes/{assembly}.fasta"
    output: "genomes/fastq/{assembly}.fq"
    message: "Using seqtk to convert {input} to FASTQ with quality score J"
    threads: 1
    shell: "seqtk seq -F 'J' {input} > {output}"

rule fragment_assemblies:
    input: "genomes/fastq/{assembly}.fq"
    output: fragsize + "genomes_fragmented/{assembly}.frag.fq.gz"
    params: config["fragment_size"]
    message: "Using seqkit to create {params}bp sliding window fragments from {input}"
    threads: 2
    shell: "seqkit sliding -j {threads} -s 1 -W {params} {input} -o {output}"

rule index_reference_bwa:
    input: ref_genome
    output: 
        idx1 = ref_genome + ".0123",
        idxamb = ref_genome + ".amb",
        idxann = ref_genome +  ".ann",
        idxbwt = ref_genome + ".bwt.2bit.64",
        idxpac = ref_genome + ".pac"
    message: "Indexing {input} with bwa-mem2"
    shell:
        """
        bwa-mem2 index {input}
        """

rule map_to_reference:
    input: 
        reference = ref_genome,
        idx1 = ref_genome + ".0123",
        idxamb = ref_genome + ".amb",
        idxann = ref_genome +  ".ann",
        idxbwt = ref_genome + ".bwt.2bit.64",
        idxpac = ref_genome + ".pac",
        query = fragsize + "genomes_fragmented/{assembly}.frag.fq.gz"
    output: temp(fragsize + "mapping/individual/{assembly}.raw.bam")
    message: "Using bwa-mem2 to map {input.query} onto {input.reference}"
    threads: 30
    params: config["bwa_parameters"]
    shell:
        """
        if [ "{threads}" = "2" ]; then
            MAPT=1
            SAMT=1
        else
            MAPT=$(awk "BEGIN {{print {threads}-int({threads}/3)}}")
            SAMT=$(awk "BEGIN {{print int({threads}/3)}}")
        fi
        bwa-mem2 mem -t $MAPT {params} -a {input.reference} {query} | samtools view -@$SAMT -F 0x04 -bS - > {output}
        """

rule sort_index_alignments:
    input: fragsize + "mapping/individual/{assembly}.raw.bam"
    output: 
        bam = fragsize + "mapping/individual/{assembly}.bam",
        idx = fragsize + "mapping/individual/{assembly}.bam.bai"
    message: "Using samtools to sort and index {input}"
    threads: 5
    shell:
        """
        samtools sort -@{threads} {input} > {output}
        samtools index {output}
        """

rule merge_alignments:
    input: expand(fragsize + "mapping/individual/{assembly}.bam", assembly = assemblies)
    output: 
        bamlist= fragsize + "mapping/individual/.bamlist",
        bam = fragsize + "mapping/alignments.bam"
    message: "Merging all of the alignments of {output.bamlist} into {output.bam}"
    threads: 30
    params: fragsize
    shell: 
        """
        ls {params}mapping/individual/*.bam > {output.bamlist}
        samtools merge -@{threads} -b {output.bamlist} -f {output} &>/dev/null
        """

rule index_alignments:
    input: fragsize + "mapping/alignments.bam"
    output: fragsize + "mapping/alignments.bam.bai"
    message: "Indexing {input} with samtools"
    threads: 1
    shell: "samtools index {input}"

rule create_popmap:
    input: fragsize + "mapping/individual/.bamlist"
    output: fragsize + "populations.map"
    message: "Creating population mapping file {output} based on genome names"
    threads: 1
    shell: 
        """
        sed 's/.bam//g' {input} | sed "s/.*\///" | paste -d' ' - <(cut -d'.' -f1 {input}) > {output}
        """

rule index_reference_samtools:
    input: ref_genome
    output: ref_genome + ".fai"
    message: "Indexing {input} with samtools"
    shell:
        """
        samtools faidx {input} > {output}
        """

rule split_regions:
    input:
        bam = fragsize + "mapping/alignments.bam",
        bai = fragsize + "mapping/alignments.bam.bai",
        fai =  ref_genome + ".fai"
    output: fragsize + "snp_discovery/reference.5kb.regions"
    message: "Splitting {input.fai} into 5kb regions for variant calling parallelization"
    threads: 1
    shell:
        """
        fasta_generate_regions.py {input.fai} 5000 > {output}
        #bamtools coverage -in {input.bam} | coverage_to_regions.py {input.fai} 500 > {output}
        """

rule call_variants:
    input:
        bam = fragsize + "mapping/alignments.bam",
        genome = ref_genome,
        regions = fragsize + "snp_discovery/reference.5kb.regions",
        populations = fragsize + "populations.map"
    output: fragsize + "variants/snps.raw.bcf"
    message: "Calling variants with freebayes"
    threads: 30
    params: config["freebayes_parameters"]
    shell: 
        """
        cat {input.regions} | parallel -k --eta -j {threads} freebayes -f {input.genome} {input.bam} \
            -C 2 --min-coverage 5 --standard-filters --populations populations.map --region {{}} \
            | vcffirstheader \
            | vcfstreamsort -w 1000 \
            | vcfuniq \
            | bcftools view -Ob - > {output}
        #freebayes-parallel {input.regions} {threads} -f {input.genome} {input.bam} -C 3 --min-coverage 5 --standard-filters {params} | bcftools view -Ob - > {output}
        """
