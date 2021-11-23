/* 
 * pipeline input parameters 
 */
params.reads = "DATA/*/*_*_{1,2}.fq.gz"

// Output location
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

// Compute and codes 
params.max_memory                 = '29.GB'
params.max_cpus                   = 4
params.max_time                   = '240.h'

// Parameter
params.adapter = '^ggcaagtgaccgtgtgtgtaaagagtgaggcgtatgaggctgtgtcggggcagaggcacaacgtttc...gcaggggagataccatgatcacgaaggtggttttcccagggcgaggcttatccattgcactccg$'
params.adapterRV = '^cggagtgcaatggataagcctcgccctgggaaaaccaccttcgtgatcatggtatctcccctgc...gaaacgttgtgcctctgccccgacacagcctcatacgcctcactctttacacacacggtcacttgcc$'
params.bwt_index = "/mnt/volume2/reference_seq"
params.overlap=6
params.gtf = "/mnt/volume2/reference_seq.gtf"
params.ref = "/mnt/volume2/reference_seq.fa"

log.info """\
         GERMS-16s - N F   P I P E L I N E    
         ===================================

         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * 	Create Channel
 */
reads_ch = Channel
			.fromFilePairs( params.reads )
			.map { 	item -> 
					sampleName = item[0];
					sampleName = sampleName.split('_')[0]
					files  = item[1];
					return [ sampleName, files ]  }
		
process fastp {
    	
	tag "$sample_id"
	label "process_medium"
	
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: '*.fq.gz'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: '*.log'
	
    input:
    tuple val(sample_id) , path(fastq) from reads_ch
	val overlap from params.overlap
	val adapter from params.adapter
	
    output:
	tuple val(sample_id) , file("${sample_id}.MERGED.fq.gz") into reads_merged_ch
	tuple val(sample_id) , file("${sample_id}.unmerged.passQC.1.fq.gz") , file("${sample_id}.unmerged.passQC.2.fq.gz") , file("${sample_id}.unmerged.failQC.1.fq.gz") , file("${sample_id}.unmerged.failQC.2.fq.gz") into reads_unmerged_ch
	file "${sample_id}.fastp.log" into log_ch_1
	
	script:
    """
	fastp --correction --merge --length_required 10 --thread $task.cpus --merged_out ${sample_id}.MERGED.fq.gz --in1 ${fastq[0]} --in2 ${fastq[1]} --out1 ${sample_id}.unmerged.passQC.1.fq.gz  --out2 ${sample_id}.unmerged.passQC.2.fq.gz --unpaired1 ${sample_id}.unmerged.failQC.1.fq.gz --unpaired2 ${sample_id}.unmerged.failQC.2.fq.gz >> ${sample_id}.fastp.log 2>&1
	echo "COMPLETED Step1 (fastp merge and correct reads pair) : ${sample_id}" >> ${sample_id}.fastp.log
	"""
}

process cutadapt {
    	
	tag "$sample_id"
	label "process_medium"
	
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: '*.fq.gz'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: '*.log'
	
    input:
    tuple val(sample_id) , file(fastq) from reads_merged_ch
	val overlap from params.overlap
	val adapter from params.adapter
	val adapterRV from params.adapterRV
	
    output:
	tuple val(sample_id) , file("${sample_id}.TRIMMED.fq.gz") into reads_trimmed_ch
	tuple val(sample_id) , file("${sample_id}.untrimmed.fq.gz") into reads_untrimmed_ch
	file "${sample_id}.cutadapt.log" into log_ch_2
	
	script:
    """
	cutadapt -j $task.cpus -a $adapter -a $adapterRV --overlap $overlap -o ${sample_id}.TRIMMED.fq.gz --untrimmed-o ${sample_id}.untrimmed.fq.gz $fastq >> ${sample_id}.cutadapt.log 2>&1
	echo "COMPLETED Step1 (cutadapt 2adapters in merged reads) : ${sample_id}" >> ${sample_id}.cutadapt.log
	"""
}

process count_reads {
    	
	tag "$sample_id"
	label "process_low"
	
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: '*.csv'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: '*.log'
	
    input:
    tuple val(sample_id) , file(fastq) from reads_trimmed_ch
	val ref from params.ref
	
    output:
	file("${sample_id}.countMatrix.csv") into sample_readCount_ch
	file "${sample_id}.counting.log" into log_ch_3
	
	script:
    """
	gunzip -c -f $fastq > reads.fq 
	python /mnt/volume2/src/counting_reads.py --input_file reads.fq --sample_id $sample_id --ref $ref >> ${sample_id}.counting.log 2>&1
	rm reads.fq
	echo "COMPLETED Step3 (Counting reads) : ${sample_id}" >> ${sample_id}.counting.log
	"""
}

