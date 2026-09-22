# Scripts to reproduce analyses in Morphological evolution and phylogeny are unlinked in Hawaiian _Doryopteris_ (Lau and Zhang et al., submitted)

## Loci analyses
### Processing raw reads into locus trees
<ol>
    <li> Align using MAFFT (web server: https://mafft.cbrc.jp/alignment/server/) </li> 
    <li> Trim gappy sites using TrimAl (`scripts/01_loci_analyses/all_samples_dory_plus_outgroup/01_trim_alignments_trimal.sh`) </li> 
    <li> Remove problematic tips (`scripts/01_loci_analyses/all_samples_dory_plus_outgroup/02_drop_problematic_tips.sh`) </li> 
    <li> Identify loci with multiple copies per sample (`scripts/01_loci_analyses/all_samples_dory_plus_outgroup/03_identify_bad_loci.R`) </li> 
    <li> Infer locus trees using IQTREE (`scripts/01_loci_analyses/all_samples_dory_plus_outgroup/04_infer_locus_trees_iqtree.sh`) </li>
    <li> Resolve cases of multiple samples per individual (manual inspection + `scripts/01_loci_analyses/all_samples_dory_plus_outgroup/05_loci_remove_sister_copy.R`) </li>
    <li> Re-align loci using MAFFT and re-infer locus trees using IQTREE </li> 
</ol> 

### Infer species tree for all _Doryopteris_ samples + outgroup
<ol>
    <li> Collapsed nodes with lower than 10 bootstrap support in each gene tree. Infer species tree using ASTRAL-II (`scripts/02_analyses/loci_run_03_astral.sh`) </li>
</ol> 

### Subset the alignments to include only Hawaiian _Doryopteris_ and _D. concolor_ samples
<ol>
    <li> Process alignments: only keep a subset of samples from HI + outgroup, retaining or dropping D. subdecipiens depending on iteration. Drop samples represented in fewer than 50% of alignments. Re-align with mafft. (`scripts/01_loci_analyses/subset_HI_plus_concolor/07_reduce_alignments.py`) </li>
    <li> Infer locus trees using IQTREE (`scripts/01_loci_analyses/subset_HI_plus_concolor/08_infer_gene_trees.sh`)</li>
    <li> Prepare ASTRAL input (`scripts/01_loci_analyses/subset_HI_plus_concolor/10_prep_astral.sh`)) </li>
    <li> Infer species tree using weighted ASTRAL and add island source to tip names (`scripts/02_analyses/12_infer_species_trees_wastral.sh` and `scripts/02_analyses/09_rename_tips_island.py`) </li>
</ol>

## SNP analyses
### Processing raw reads into SNP data
<ol>
    <li> tbc </li>
</ol>

### Infer population-level relationships
<ol>
<li> Perform PCA using SNPRelate in R
    <li> Pacific _Doryopteris_ excluded (`scripts/02_analyses/snp_run_01_SNPRelate_HI.R`) </li>
    <li> Pacific _Doryopteris_ included (`scripts/02_analyses/snp_run_02_SNPRelate_HI_plus_outgroup.R`) </li>
</li>
<li> Visualized phylogenetic network using SplitsTree4 (run using GUI; https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) </li>
</ol>

## Plot figures
<ol>
    <li> Figure 1: Distribution of _Doryopteris_ (`scripts/03_plotting/plot_fig1_map_ranges.R`) </li>
    <li> Figure 2: Backbone phylogeny (`scripts/03_plotting/plot_fig2_astral.R`) </li>
    <li> Figure 3: Phylogenies of Hawaiian _Doryopteris_ (`scripts/03_plotting/plot_fig3_wastral.R`) </li>
    <li> Figure 4: PCA and phylogenetic network (results of PCA visualized using `scripts/03_plotting/plot_fig4_snprelate.R`; results of SplitsTree4 directly available from the software) </li>
</ol>

Adobe Illustrator and Inkscape were used to modify the aesthetic aspect of the figures.

## Supplmentary materials
<ol>
    <li> Infer MatK phylogeny using IQTREE (`scripts/04_suppmat/matK_run_iqtree.sh`) </li>
</ol>

## Metadata
<ol>
    <li> `metadata/combined_sample_data.csv` </li>
</ol>
