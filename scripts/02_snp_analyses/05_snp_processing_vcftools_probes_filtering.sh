

##############################
### final filtering parameters
##############################
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --gzvcf Doryopteris_flanking_SNPs.vcf.gz --min-alleles 2 --max-alleles 2 --max-missing 0.5 --minDP 3 --maxDP 100 --recode --recode-INFO-all --out Doryopteris_initial_filtering_all2
#vcftools --gzvcf Doryopteris_probes_SNPs_variant_only.vcf.gz --max-missing 0.7 --min-meanDP 3 --max-meanDP 100 --recode --recode-INFO-all --out Doryopteris_initial_filtering
#After filtering, kept 53999 out of a possible 391059 Sites

#gives missing proportion of loci for each individual
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf Doryopteris_initial_filtering_all.recode.vcf --missing-indv

#create a list of individuals with at least 50% missing data
#awk '$5 > 0.5' out.imiss | cut -f1 > lowDP50.indv


#gives missing proportion of loci for each individual

#keep samples for two different data sets
#HI only
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf Doryopteris_initial_filtering_all.recode.vcf --keep HI_only.txt --min-alleles 2 --max-alleles 2 --max-missing 0.5 --minDP 3 --maxDP 100 --recode --recode-INFO-all --out Ready_for_plink_HI

/programs/plink-1.9-x86_64-beta7/plink --vcf Ready_for_plink_HI.recode.vcf --make-bed --allow-extra-chr --double-id --out pop_sorted --set-missing-var-ids @:#

#tests for pairwise LD in a sliding window (5kb window with a step size of 0.5 and r^2 of 0.5); this value will change based on the plateau from PopLDdecay
/programs/plink-1.9-x86_64-beta7/plink --bfile pop_sorted --indep-pairwise 20 10 0.2 --allow-extra-chr --double-id --threads 4

#this grabs the correct positions, but chr1-20 need to be modified in the To_prune.txt to put the prefix in
awk '{gsub(":", " ");print}' plink.prune.out > To_prune.txt

#/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools to prune out selected SNPs from plink
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf Ready_for_plink_HI.recode.vcf --exclude-positions To_prune.txt --recode --recode-INFO-all --out HI_only_LD_pruned_SNPs

#gives missing proportion of loci for each individual
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf HI_only_LD_pruned_SNPs.recode.vcf --missing-indv

#average depth for each individual
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf HI_only_LD_pruned_SNPs.recode.vcf --depth 

#observed and expected heterozygosity
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf HI_only_LD_pruned_SNPs.recode.vcf --het


#create fasta file for downstream analyses
/local/workdir/jbl256/Installed_programs/vcf2phylip-2.8/vcf2phylip.py --input HI_only_LD_pruned_SNPs.recode.vcf -f -n -b

mv Ready_for_plink_HI.recode.vcf HI_only_linked_SNPs.recode.vcf 
/local/workdir/jbl256/Installed_programs/vcf2phylip-2.8/vcf2phylip.py --input HI_only_linked_SNPs.recode.vcf -f -n -b


#######
# HI plus outgroup
#######

#HI+outgroups only
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf Doryopteris_initial_filtering_all.recode.vcf --keep HI_outgroup.txt --min-alleles 2 --max-alleles 2 --max-missing 0.5 --minDP 3 --maxDP 100 --recode --recode-INFO-all --out Ready_for_plink_HI_plus_out

/programs/plink-1.9-x86_64-beta7/plink --vcf Ready_for_plink_HI_plus_out.recode.vcf --make-bed --allow-extra-chr --double-id --out pop_sorted --set-missing-var-ids @:#

#tests for pairwise LD in a sliding window (5kb window with a step size of 0.5 and r^2 of 0.5); this value will change based on the plateau from PopLDdecay
/programs/plink-1.9-x86_64-beta7/plink --bfile pop_sorted --indep-pairwise 20 10 0.2 --allow-extra-chr --double-id --threads 4

#this grabs the correct positions, but chr1-20 need to be modified in the To_prune.txt to put the prefix in
awk '{gsub(":", " ");print}' plink.prune.out > To_prune.txt

#/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools to prune out selected SNPs from plink
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf Ready_for_plink_HI_plus_out.recode.vcf --exclude-positions To_prune.txt --recode --recode-INFO-all --out HI_plus_outgroups_LD_pruned_SNPs

#gives missing proportion of loci for each individual
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf HI_plus_outgroups_LD_pruned_SNPs.recode.vcf --missing-indv

#average depth for each individual
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf HI_plus_outgroups_LD_pruned_SNPs.recode.vcf --depth 

#observed and expected heterozygosity
/local/workdir/jbl256/Installed_programs/vcftools-0.1.16/bin/vcftools --vcf HI_plus_outgroups_LD_pruned_SNPs.recode.vcf --het


#create fasta file for downstream analyses
/local/workdir/jbl256/Installed_programs/vcf2phylip-2.8/vcf2phylip.py --input HI_plus_outgroups_LD_pruned_SNPs.recode.vcf -f -n -b

mv Ready_for_plink_HI_plus_out.recode.vcf HI_plus_outgroups_linked_SNPs.recode.vcf 
/local/workdir/jbl256/Installed_programs/vcf2phylip-2.8/vcf2phylip.py --input HI_plus_outgroups_linked_SNPs.recode.vcf -f -n -b
