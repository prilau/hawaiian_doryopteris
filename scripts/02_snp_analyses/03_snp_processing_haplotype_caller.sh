mkdir Duplicates_marked
mkdir GVCF
mkdir Completed_BAMs_after_SNPs

#create sequence dictionary in Picard (needed for GATK analysis), need to already have the genome indexed with BWA

for file in sorted_bam_files/*.bam
do
echo "Marking duplicates in $file "
name=`basename $file .bam`
/home/jbl256/Installed_programs/gatk-4.3.0.0/gatk MarkDuplicates -I $file -O Duplicates_marked/$name.duplicates.bam -M Duplicates_marked/$name.dup_metrics.txt
done


#index bam files before calling SNPs
for file in Duplicates_marked/*.bam
do
echo "Indexing $file "
name=`basename $file .bam`
/home/jbl256/Installed_programs/samtools-1.16.1/samtools index $file
done

#Use HaplotypeCaller for each sample
for file in Duplicates_marked/*.bam
do
	echo "Calling SNPs on $file "
	name=`basename $file .bam`
	/home/jbl256/Installed_programs/gatk-4.3.0.0/gatk HaplotypeCaller -R ../reference_files/full_locus_only_one_target_accession.fasta -I $file -O GVCF/$name.g.vcf.gz -ERC GVCF
	mv $file Completed_BAMs_after_SNPs
	mv $file.bai Completed_BAMs_after_SNPs
done



