# Phylogenetic tree inference of Hawaiian Doryopteris
This repository contains scripts to reproduce the analyses in the manuscript 'Morphological evolution and phylogeny are unlinked in Hawaiian _Doryopteris_' (Lau and Zhang et al., submitted). Specifically, our analyses are based on two types of data we generated from raw reads - (i) alignments of 316 loci (exon regions) and (ii) single nucleotide polymorphic (SNP) sites. 
Raw sequences are deposited at xxx. Alignments and the SNP dataset can be found in Zenodo repository xxx.

### Analyses based on alignments of loci
In the first part (`scripts/01_loci_analyses/all_samples_dory_plus_outgroup`), we processed the raw reads into alignments and used ASTRAL-II to infer a species phylogeny that includes all Doryopteris samples. In the second part (`scripts/01_loci_analyses/subset_HI_plus_concolor`), we subset the alignments to keep only Hawaiian samples and their sister group. We ran two weighted ASTRAL analyses, one excluding *D. subdecipiens* and the other including those samples.

### SNP analyses
Next, we produced a SNP dataset to infer the population-level relationships (`scripts/02_snp_analyses`). We performed two principal component analyses using the R package SNPRelate (Zheng et al., 2012; https://github.com/zhengxwen/SNPRelate), one with Hawaiian samples only, and the other one with both Hawaiian samples and Pacific *D. concolor*. We also used SplitTree4 (Hudson & Bryant 2006; https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) to visualize phylogenetic networks (run using GUI; ).

### Plot figures
We plotted most parts of the figures in R (see `scripts/03_plotting`), and used Adobe Illustrator and Inkscape to modify the aesthetic aspect.


* Figure 1: Distribution of _Doryopteris_ (`plot_fig1_map_ranges.R`)
* Figure 2: Backbone phylogeny (`plot_fig2_astral.R`)
* Figure 3: Phylogenies of Hawaiian _Doryopteris_ (`plot_fig3_wastral.R`)
* Figure 4: PCA and phylogenetic network (results of PCA visualized using `plot_fig4_snprelate.R`; results of SplitsTree4 directly available from the software)



## Supplmentary materials
We also inferred a tree using the plastid gene matK. (`scripts/04_suppmat/matK_run_iqtree.sh`)
