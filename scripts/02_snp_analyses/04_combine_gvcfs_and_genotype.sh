

#combine individual gVCFs into one file and call SNPs on combined gVCF file
/home/jbl256/Installed_programs/gatk-4.3.0.0/gatk GenomicsDBImport -R ../../reference_files/full_locus_only_one_target_accession.fasta -V samples.list -L interval.list --genomicsdb-workspace-path flanking_regions --batch-size 20
/home/jbl256/Installed_programs/gatk-4.3.0.0/gatk GenotypeGVCFs -R ../../reference_files/full_locus_only_one_target_accession.fasta -V gendb://flanking_regions --output Doryopteris_flanking_SNPs.vcf.gz --include-non-variant-sites TRUE

/home/jbl256/Installed_programs/gatk-4.3.0.0/gatk GenotypeGVCFs -R ../../reference_files/full_locus_only_one_target_accession.fasta -V gendb://flanking_regions --output Doryopteris_flanking_SNPs_variant_only.vcf.gz --include-non-variant-sites FALSE
