# Scripts to reproduce analyses in Morphological evolution and phylogeny are unlinked in Hawaiian _Doryopteris_ (Lau and Zhang et al., submitted)

### Processing raw reads into locus trees
<ol>
<li> Align using MAFFT (web server: https://mafft.cbrc.jp/alignment/server/) </li> 
<li> Trim gappy sites using TrimAl (`scripts/01_data_processing/loci_process_01_trimal.sh`) </li> 
<li> Identify good and bad loci (`scripts/01_data_processing/loci_process_02_good_loci.sh`) </li> 
<li> Search for best-fit substitution models (`scripts/02_analyses/loci_run_01_modelfinder.sh`) </li>
<li> Infer locus trees using IQTREE2 (`scripts/02_analyses/loci_run_02_iqtree.sh`) </li>
<li> Resolve cases of multiple samples per individual (manual inspection) </li>
<li> Infer locus trees using IQTREE2 (`scripts/02_analyses/loci_run_02_iqtree.sh`) </li>
</ol> 

### 2. Infer species tree
<ol>
<li> Infer species tree using all _Doryopteris_ and outgroup samples using ASTRAL-II (`scripts/02_analyses/loci_run_03_astral.sh`) </li>
</li> Infer species tree using all Hawaiian _Doryopteris_ samples and Pacific _Doryopteris_ samples using weighted ASTRAL (`scripts/02_analyses/loci_run_04_wastral.sh`)
    </li> _Doryopteris subdecipiens_ samples excluded </li>
    </li> _Doryopteris subdecipiens_ samples included </li>
</li>
</ol>

### Processing raw reads into SNP data
i. tbc

### Infer population-level relationships
<ol>
<li> Perform PCA using SNPRelate in R
    <li> Pacific _Doryopteris_ excluded (`scripts/02_analyses/snp_run_01_SNPRelate_HI.R`) </li>
    <li> Pacific _Doryopteris_ included (`scripts/02_analyses/snp_run_02_SNPRelate_HI_plus_outgroup.R`) </li>
</li>
<li> Visualized phylogenetic network using SplitsTree4 (run using GUI; https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) </li>
</ol>

### Plot figures
<ol>
<li> Figure 1: Distribution of _Doryopteris_ (tbc) </li>
<li> Backbone phylogeny (`scripts/03_plotting/loci_plot_astral.R`) </li>
<li> Phylogenies of Hawaiian _Doryopteris_ (tbc) </li>
<li> PCA and phylogenetic network (results of PCA visualized using `scripts/03_plotting/snp_plot_snprelate.R`; results of SplitsTree4 directly available from the software) </li>
</ol>

Adobe Illustrator and Inkscape were used to modify the aesthetic aspect of the figures.
