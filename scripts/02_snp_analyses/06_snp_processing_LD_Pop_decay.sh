#############################################
##### LD pop decay    ########
#############################################
#calculate LD on all samples
#/home/jbl256/Installed_programs/PopLDdecay/bin/PopLDdecay -InVCF Doryopteris_initial_filtering_all.recode.vcf -OutStat LDdecay_Dory -MaxDist 10 -MAF 0.05 -OutType 1
/home/jbl256/Installed_programs/PopLDdecay/bin/PopLDdecay -InVCF Doryopteris_initial_filtering_all.recode.vcf -OutStat LDdecay_Dory -MaxDist 5 -OutType 1

#plot figure for Rice genotypes
perl /home/jbl256/Installed_programs/PopLDdecay/bin/Plot_OnePop.pl -inFile LDdecay_Dory.stat.gz -output Fig_LDdecay_Dory_all_samples -keepR



