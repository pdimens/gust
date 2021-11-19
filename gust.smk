import glob
configfile: "config.yaml"
ref_genome = config["reference_genome"]

rule all:
    input: "variants/snps.raw.bcf"


rule fastq_convert:
    input:  "genomes/{assembly}.fasta"
    output: "genomes/fastq/{assembly}.fq"
    message: "Using seqtk to convert {input} to FASTQ format with dummy quality score J"
    threads: 1
    shell: "seqtk seq -F 'J' {input} > {output}"

rule fragment_assemblies:
    input: "genomes/fastq/{assembly}.fq"
    output: "genomes_fragmented/{assembly}.frag.fq.gz"
    message: "Using seqkit to create 50bp sliding window fragments on {input}"
    threads: 2
    shell: "seqkit sliding -j {threads} -s 1 -W 50 {input} -o {output}"

rule index_reference:
    input: ref_genome
    output: ref_genome + "fai"
    message: "Indexing {input}"
    shell: "bwa-mem2 index {input}"

rule map_to_reference:
    input: 
        reference = ref_genome,
        refindex = ref_genome + ".fai",
        query = "genomes_fragmented/{assembly}.frag.fq.gz"
    output: temp("mapping/individual/{assembly}.raw.bam")
    message: "Using bwa-mem2 to map {input.query} onto {input.reference}"
    threads: 20
    params:
        mapthreads = {threads} - ({threads} // 3),
        samthreads = {threads} // 3
    shell:
        """
        bwa-mem2 mem -t {params.mapthreads} -a {input.reference} {query} | samtools view -@{params.samthreads} -F 0x04 -bS - > {output}
        """

rule process_bamfiles:
    input: "mapping/{assembly}.raw.bam"
    output: 
        bam = "mapping/individual/{assembly}.bam",
        idx = "mapping/individual/{assembly}.bam.bai"
    message: "Using samtools to sort and index {input}"
    threads: 5
    shell:
        """
        samtools sort -@{threads} {input} > {output}
        samtools index {output}
        """

rule merge_alignments:
    input: glob.glob('mapping/individual/*.bam')
    output: 
        bamlist= "mapping/individual/.bamlist",
        bam = "mapping/alignments.bam"
    message: "Merging all of the alignments of {output.bamlist} into {output.bam}"
    threads: 20
    shell: 
        """
        ls mapping/individual/*.bam > {output.bamlist}
        samtools merge -@{threads} -b {output.bamlist} -f {output} &>/dev/null
        """

rule index_alignments:
    input: "mapping/alignments.bam"
    output: "mapping/alignments.bam.bai"
    message: "Indexing {input}"
    threads: 1
    shell: "samtools index {input}"

rule create_popmap:
    input: "mapping/individual/.bamlist"
    output: "misc/populations.map"
    message: "Creating population mapping file {output} based on genome names"
    threads: 1
    shell: 
        """
        sed 's/.bam//g' {input} | sed "s/.*\///" | paste -d' ' - <(cut -d'.' -f1 {input})
        """

rule assess_coverage:
    input:
        bam = "mapping/alignments.bam",
        bai = "mapping/alignments.bai",
        fai =  ref_genome + ".fai"
    output: "snp_discovery/reference.5kb.regions"
    message: "Splitting {input.fai} into 5kb regions for variant calling parallelization"
    threads: 1
    shell:
        """
        fasta_generate_regions.py {input.fai} 5000 > {output}
        #bamtools coverage -in {input.bam} | coverage_to_regions.py {input.fai} 500 > {output}
        """

rule call_variants:
    input:
        bam = "mapping/alignments.bam",
        genome = ref_genome,
        fai = ref_genome + ".fai",
        regions = "snp_discovery/reference.5kb.regions",
        populations = "misc/populations.map"
    output: "variants/snps.raw.bcf"
    message: "Calling variants with freebayes"
    threads: 20
    shell: 
        """
        freebayes-parallel {input.regions} {threads} -f {input.genome} {input.bam} -C 3 --min-coverage 5 --standard-filters | bcftools view -Ob - > {output}
        """
