mkdir SAM_files
mkdir bam_files
mkdir sorted_bam_files

#to index a fasta file
/home/jbl256/Installed_programs/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index full_locus_only_one_target_accession.fasta

#Read group information starts with "@RG
#ID: is unique identifier of the samples, for now doing the sample name and the barcode info
#SM: is the sample name
#PL: is the sequencing equipment, in almost all cases this will be Illumina
#PU: is the run identifier, the lane, followed by the specific barcode of the sample
#LB: is the library count

for file in /data/jbl256/Doryopteris/cleaned_reads/*R1.fastq.gz
do
	name=`basename $file .R1.fastq.gz`
	echo "Mapping reads for $name to reference"
	forward=$name".R1.fastq.gz"
	reverse=$name".R2.fastq.gz"
	design=`basename $file .R1.fastq.gz`
	#perform the alignment
	/home/jbl256/Installed_programs/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 8 -R "@RG\tID:$design.run1\tSM:$design\tPL:IlluminaNovaSeq\tPU:HTNMKDSXX\tLB:Dory_flanking" ../reference_files/full_locus_only_one_target_accession.fasta ../cleaned_reads/$forward ../cleaned_reads/$reverse > SAM_files/$name.sam
done

#convert SAM to BAM for sorting

for file in SAM_files/*.sam
do
	echo "Convert $file to to BAM"
	name=`basename $file .sam`
	/home/jbl256/Installed_programs/samtools-1.16.1/samtools view -S -b $file > bam_files/$name.bam
	rm $file
done

#Sort BAM for SNP calling
for file in bam_files/*.bam
do
	echo "Sort $file"
	name=`basename $file .bam`
	readid=$name
	/home/jbl256/Installed_programs/samtools-1.16.1/samtools sort -o sorted_bam_files/$readid.bam $file
	rm $file
done

