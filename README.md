# pheromoneFinder
<i>Written by Sri Srikant, as part of graduate research in the labs of Andrew Murray and Rachelle Gaudet.</i>

Code-base to identify fungal pheromones from genome sequences. Written for publication being prepared from the Murray Lab. Started writing the code to identify the pheromone of Y. lipolytica to inform my project on the pheromone exporting ABC exporter in Y. lipolytica.

Code:
1. scripts/pheromoneFinder_*yyyymmdd*.py:	Python scripts written to be used with SLURM scheduler on Odyssey cluster at Harvard. Functions within the code provide the methods to be run in sequence to identify pheromone candidates from fungal genome sequences. Primarily used with genomes of yeasts in the Saccharomycotina clade of Ascomyctoa.
2. notebooks/matingGenePromoter_*yymmdd*.ipynb:	Python notebook written in python 2.7 to test code fragments that are incorporated into pheromoneFinder.py scripts. Other code-cells are written to analyze outputs from third-party binaries and plot results.
3. notebooks/pheromoneFinder_py*2/3*-testbed.ipynb:	Jupyter notebooks that are for testing code-fragments. python 3 notebook is used to identify clades within specific time-horizons using the phylogenetic tree from Shen XX et al, Cell 2018. 
