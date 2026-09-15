# Scripts to reproduce analyses in Morphological evolution and phylogeny are unlinked in Hawaiian _Doryopteris_ (Lau and Zhang et al., submitted)

### Processing raw reads into locus trees
<ol>
<li> Align using MAFFT and trim gappy sites using TrimAl (`scripts/loci_process_muscle_trimal.sh`) </li> 
<li> Search for best-fit substitution models (`scripts/loci_run_modelfinder.sh`) </li>
<li> Infer locus trees using IQTREE2 (`scripts/loci_run_iqtree.sh`) </li>
<li> Resolve cases of multiple samples per individual (manual inspection) </li>
<li> Infer locus trees using IQTREE2 (`scripts/loci_run_iqtree.sh`) </li>
</ol> 

### 2. Infer species tree
<ol>
<li> Infer species tree using all _Doryopteris_ and outgroup samples using ASTRAL-II (`scripts/loci_run_astral.sh`) </li>
</li> Infer species tree using all Hawaiian _Doryopteris_ samples and Pacific _Doryopteris_ samples using weighted ASTRAL (`scripts/loci_run_wastral.sh`)
    </li> _Doryopteris subdecipiens_ samples excluded </li>
    </li> _Doryopteris subdecipiens_ samples included </li>
</li>
</ol>

### Processing raw reads into SNP data
i. tbc

### Infer population-level relationships
<ol>
<li> Perform PCA using SNPRelate in R
    <li> Pacific _Doryopteris_ excluded (`scripts/loci_run_SNPRelate_HI.sh`) </li>
    <li> Pacific _Doryopteris_ included (`scripts/loci_run_SNPRelate_HI_plus_outgroup.sh`) </li>
</li>
<li> Visualized phylogenetic network using SplitsTree4 (run using GUI; https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) </li>
</ol>

### Plot figures
<ol>
<li> Figure 1: Distribution of _Doryopteris_ (distribution retrieved using `scripts/plot_map.sh`) </li>
<li> Backbone phylogeny (`scripts/plot_astral.sh`) </li>
<li> Phylogenies of Hawaiian _Doryopteris_ (`scripts/plot_wastral.sh`) </li>
<li> PCA and phylogenetic network (results of PCA visualized using `scripts/plot_SNPRelate.sh`; results of SplitsTree4 directly available from the software) </li>
</ol>

Adobe Illustrator and Inkscape were used to modify the aesthetic aspect of the figures.
