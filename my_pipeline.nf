nextflow.enable.dsl = 2

params.outdir = "Next_ik"

process FASTQC{
	tag "FASTQC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${params.outdir}/QC", pattern: 'fastqc_out/*.html', mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
	path "fastqc_out/*"

	script:
	""" 
	mkdir fastqc_out
	fastqc $reads -o fastqc_out
	"""
}

process CUTADAPT{
	tag "CUTADAPT on $name using $task.cpus CPUs and $task.memory memory" 
	publishDir "${params.outdir}/trimmed/", mode:'copy'
		
	input:
	tuple val(name), path(reads)
	
	output:
	tuple val("trimmed_${name}"), path("trimmed_*.fastq.gz")

	script:
	"""
	mkdir trimmed
	cutadapt -m 10 -q 20 -j 8 -o trimmed_${reads[0]} -p trimmed_${reads[1]} ${reads}
	"""
}

process BWAMEM {
	tag "BWAMEM on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/mapped", mode:'copy'

	input:
	tuple val(name),path(reads)
	
	output:
	tuple val(name), path("*.sam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""

	"""
	bwa mem -R ${rg} -t $task.cpus -o ${name}.sam ${params.refindex} $reads  
	"""
}

process SAMCLN {
	tag "SAMCLN on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/bams/", mode:'copy'

	input:
	tuple val(name),path(reads)
	
	output:
	
	path "*markdup*"

	script:
	"""
	picard SortSam INPUT=$reads OUTPUT=${name}.srt.bam SORT_ORDER=coordinate
	picard MarkDuplicates INPUT=${name}.srt.bam OUTPUT=${name}.srt.markdup.bam METRICS_FILE=${name}.markdup.txt
	picard BuildBamIndex INPUT=${name}.srt.markdup.bam
	"""
}

process VARCAL {
	tag "VARCAL on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/VCFs/", mode:'copy'

	input:
	path bam_files
	
	output:
	
	path "${bam_files.getSimpleName()}.vcf*"

	script:
	"""
	gatk HaplotypeCaller -R ${params.refvcf} -I ${bam_files} -O ${bam_files.getSimpleName()}.vcf
	"""
}

rawfastq = channel.fromFilePairs("${params.outdir}/raw_fastq/wVA*_R{1,2}*", checkIfExists: true)
//rawfastq_ch = channel.fromPath("${params.outdir}/raw_fastq/w*", checkIfExists: true) //creates tuple with [name,path/name]
//		.map{[it.getName(),it]}

workflow {
	FASTQC(rawfastq)
	trimmed = CUTADAPT(rawfastq)
	sams = BWAMEM(trimmed)
	bams = SAMCLN(sams)
bams_filt = bams
	.flatten()
	.filter(~/.*.bam/)    
		
	VARCAL(bams_filt)
}
