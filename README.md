```
    _______.  ______  _______   __   _______  _______ ___   
    /       | /      ||       \ |  | |   ____||   ____|__ \  
   |   (----`|  ,----'|  .--.  ||  | |  |__   |  |__     ) | 
    \   \    |  |     |  |  |  ||  | |   __|  |   __|   / /  
.----)   |   |  `----.|  '--'  ||  | |  |     |  |     / /_  
|_______/     \______||_______/ |__| |__|     |__|    |____| 

```
                                                           
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Table of Contents
1. [INTRODUCTION](#INTRODUCTION)
2. [PREREQUISITES](#PREREQUISITES)
3. [INSTALLATION](#INSTALLATION)
4. [USAGE](#USAGE)
5. [EXAMPLE](#EXAMPLE)
6. [GUIDELINE](#GUIDELINE)
7. [UTILS](#UTILS)
8. [CREDITS](#CREDITS)
9. [LICENSE](#LICENSE)
10. [CONTACT](#CONTACT)


# INTRODUCTION 
This is the next version of our scdiff software suite (https://github.com/phoenixding/scdiff). 
scdiff was proven successful in inferring cell differentiation trajectories and the underlying regulatory networks 
(e.g,https://doi.org/10.1016/j.stem.2018.09.009 and https://doi.org/10.1016/j.stem.2019.12.009). 
However, with the rapid development of single-cell technologies, many new computational challenges arised 
in reconstructing the cell dynamics (dynamic trajectories and gene regulation) in various biological processes. 
To address those challenges, we would require more efficient, scalable, accurate and interactive models for better
exploring the expensive single-cell genomics data. scdiff2 is developed to meet the emergent needs in the field for such methods.  



## 1. New scdiff2 now handle huge single cell data efficiently! 
As the scale of the single-cell RNA-seq datasets is ever-increasing (from hundred cells ->tens of thousand cells and even more. 
The memory and time efficiency for the original scdiff is becoming a bottleneck of its application.  
Here, we have been developing the next version, scdiff2, that uses HDF5, Sparse-matrix, and multi-threading techniques to reduce the resource requirement of the program while improving the efficiency.  Besides, we also incorporated many popular clustering and trajectory methods (mostly implemented by scanpy https://scanpy.readthedocs.io/en/stable)
in a "prerun" program to learn an initial trajectory for the future PGM refinement (by HMM-like Probabilistic graphical models).    
scdiff2 now can finish processing 40k cells (~10k genes/cell) within 1 hour @ a desktop: Ryzen 3500 6 cores, 16G RAM) with PGM refinement (--ncores 10 --llhCut 0.05).
Without the PGM refinement, it can complete in a few minutes (e.g, 7 mins using --ncores 10 --maxloop 10 parameters).   


## 2. New scdiff2 now is fully customizable!  
First, the selection of the root node (cells) is critical for the tree-structure cell trajectory inference. 
In scdiff2, we combined the trajectory from PAGA (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x)
and the capture time (the actual time of the cells) to infer a potential root of the tree.  
Based on the inferred root, we will build the trajectory based on both gene expression and cell capture time information.
The users are also allowed to customize the root. Any clusters are allowed to be set as the root, the program will automatically 
learn the best tree structure (and underlying regulatory networks) accordingly.     

Second, the new scdiff2 was composed of many moving pieces now, each can be customized and modified (via the provided Graph class in scdiff2, will provide the tutorial soon).  
It has two passes: in the first pass, we learn an initialization using a provided prerun program (any user-provided);
in the second pass, we used a Probablistic Graphical Model (PGM) based strategy (described in the original scdiff paper) 
to iteratively refine the trajectory (tree) from the prerun initialization. 

## 3 New scdiff2 inherits the interactive visualization functionality!
Such interactive capability was proven to be very useful in finding key regulatory factors (TFs), genes (DE genes), and functions (GO term analysis) in our previous
collaborations with developmental biologists and cancer biologists. We are committed to keep improving the interactive functionality as per requests from the users.  


# PREREQUISITES
* Python3.6+
It was installed by default for most Linux distribution and MAC.  
If not, please check [https://www.python.org/downloads/](https://www.python.org/downloads/) for installation 
instructions.   
__PLEASE NOTE__: python 2.7 is no longer supported by the scdiff software suite. Please consider upgrading to Python3.6+  

* Python packages dependencies:    
	-- scipy>0.13.3    
	-- numpy>1.8.2   
	-- scikit-learn>=0.20,<=0.22  
	-- matplotlib>=3.1.2  
	-- imbalanced_learn<0.5.0  
	-- anndata>=0.7  
	-- scanpy>=1.5  
	-- pandas>=0.23  
	-- h5py>=2.10  
The python setup.py script (or pip) will try to install these packages automatically.
Also, the pip3 installation will install those libraries automatically.  However, please install them manually if, by any reason, the automatic 
installation fails. 


# INSTALLATION
 There are 2 options to install scdiff.  
* __Option 1: Install from download directory__   
	cd to the downloaded scdiff package root directory

	```shell
	$ cd scdiff
	```
	run python setup to install   

	```shell
	$ python3 setup.py install
	```
	
* __Option 2: Install from Github__ (recommended):    

	python 3: 
	```shell
	$ sudo pip3 install --upgrade https://github.com/phoenixding/scdiff2/zipball/master
	```
		
# USAGE
(1) 1st pass: prerun  
```shell
usage: prerun [-h] -i INPUT -o OUTPUT -f FORMAT

scdiff2 pre-run

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input single cell RNA-seq expression data
  -o OUTPUT, --output OUTPUT
                        output directory
  -f FORMAT, --format FORMAT
                        the format of input expression, either raw/norm (raw:
                        raw read counts, norm: normalized expression
                        
```
For the prerun, it has 3 parameters as the following:     
	 -- -i expression_matrix  
		expression_matrix is the same as the orignal scdiff (https://github.com/phoenixding/scdiff#inputs-and-pre-processing).  
		1st column: cell id  
		2nd column: capture time (must real value)    
		3nd column: label (label for each cell in visualization, not used in the inference, can be set as 'NA' if unknown)    
		4th- columns: expression for all genes  
	-- -o outputdir   
	outputdir specifies the output directory.  
	-- -f format  
	It can be set as 'raw' (represent raw reads expression) or 'norm' (represents log normalized expression)
	'raw' format should be good for all inputs (best for the raw reads expression matrix, but should be also fine for normalized expression).  
	If setting the 'norm', no more processing (e.g, filtering cells, genes and normalization, batch effect remover) will be performed.   
	The program will directly use the supplied expression matrix as it is, which is not good in many cases, especially when you haven't normalized the data properly. 
	

(2) 2nd pass: scdiff2
```
usage: scdiff2.py [-h] -i INPUT -o OUTPUT -t TFDNA [--etfListFile ETFLISTFILE]
                  [--log2fc LOG2FC] [--ncores NCORES] [--root ROOT]
                  [--llhCut LLHCUT] [--maxloop MAXLOOP]

scdiff2 main

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        h5ad result from pre-run
  -o OUTPUT, --output OUTPUT
                        output directory
  -t TFDNA, --tfdna TFDNA
                        TF-DNA interaction data
  --etfListFile ETFLISTFILE
                        By default, this program recognizes 1.6k TFs
                        (collected in human and mouse). Users are able to
                        provide a customized list of TFs　using this option
                        (e.g, for another species).
  --log2fc LOG2FC       By default, scdiff uses log2 Fold change
                        0.6(=>2^0.6~=1.5) as the cutoff for differential genes
                        (together with student t-test p-value cutoff 0.05).
                        Users can customize this log fold change cutoff.
  --ncores NCORES       # of allocated cpu cores for multi-threading the job
                        (4 by default)
  --root ROOT           Set the root (of the tree) as an input cluster ID
                        (e.g., 0 from the prerun result)
  --llhCut LLHCUT       The convergence likelihood cutoff, the program stops
                        if the cell assignment likelihood improvement is
                        smaller than this cutoff (e.g. 0.05: 5 percent)
  --maxloop MAXLOOP     The max # of loops allowed for the PGM based iterative
                        refinment. Set it to 0 to directly use the clustering
                        and trajectory results from the prerun program (scanpy
                        based). Only the regulatory networks (TFs) and the
                        interactive visulzation page will be learned and
                        generated

```
For the scdiff2 main program, it has 3 *REQUIRED* parameters:   
	-- -i h5ad result  
		the results from the first run  
	-- -o outputdir
		outputdir specifies the output directory.  
	-- -t tf_dna
		tf_dna specifies the tf-dna interaction file, it's the same format as the original scdiff.   
	all other parameters are optional.   
	
an example h5ad looks like the following:
```
AnnData object with n_obs × n_vars = 152 × 5342
    obs: 'label', 'time', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'leiden'
    var: 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'diffmap_evals', 'leiden', 'leiden_colors', 'leiden_sizes', 'neighbors', 'paga', 'pca', 'rank_genes_groups', 'umap'
    obsm: 'X_diffmap', 'X_pca', 'X_umap'
    varm: 'PCs'
    obsp: 'connectivities', 'distances'
```
The scdiff2 main programs used the following attributes, make sure to include them if you want to supply scdiff2 with your customized h5ad result. 
```
 AnnData object with n_obs × n_vars = 152 × 5342
    obs: 'time', 'leiden'
    var: 
    uns:  'paga'
    obsm: 'X_diffmap', 'X_pca', 'X_umap'
``` 

The --input for the scdiff2 must be in anndata format(https://anndata.readthedocs.io/en/stable/anndata.AnnData.html). 
You can use the provided 'prerun' program or 'scanpy' for the prerun analysis.

# EXAMPLE
## (1) Example input files   
We have provided an example (inputs and outputs) in the example directory [example/](example)  

a) example expression file for pre-run : [example/example.E](example/example.E)  
b) example　tf_dna file: [example/Human_TF_targets.txt](example/Human_TF_targets.txt)  

After installing scdiff2 successfully, cd to the example folder and run the following commands to analyze the example single-cell genomics data. 
```shell
$ prerun -i example.E -o example_out -f raw
$ scdiff2 -i example_out/example.E.h5ad -o example_out -t example_tfdna.txt --ncores 10
```

## (2) Example output files    
You can check all the output results under the output directory:       
i) Result with the PGM refinement (--llhCut -0.05 --maxloop 5) : [example/example_out](example/example_out)  
ii Result *without* the PGM refinment (--maxloop 0) : [example/example_out_no_PGM_refinement](example/example_out_no_PGM_refinement)   

a) A PGM refined tree-structure cell differentiation (MEF-> Neuron) trajectory for the example data (from Treutlein et al. 2016 Nature Paper on Neuron reprogramming).    
PGM refined clustering (colored by scdiff2 clusters and colors)   
![example/example_out/figures/umap_clustering.jpg](example/example_out/figures/umap_clustering.jpg)  
b) A PGM refined PAGA connectivity tree  
![iamges/paga_Traj.jpg](images/paga_Traj.jpg)  
c) A PGM refined scdiff tree model with the regulatory information (TFs) (interactive visualized).   
![images/scdiff2_model.jpg](images/scdiff2_model.png)
 

PGM VS. no refinement performance comparison (for this example data): 

__Likelihood__:    
With PGM refinement (-1196.0434282645167)> Without PGM refinement (-1334.8432584323705) => Total improvement 10.34%  

__Running time__:    
With PGM refinement (9mins) > without PGM refinement (1min) => 9 times slower   

# GUIDELINE  
### Hardware    
For best experience, please use >16G RAM and try to use any many threads (--ncores) as possible.  
Increasing --ncores (# of allocated cpu cores) may slightly increase the memory cost (limit --ncores if the memory resource is limited).
An a guideline, a 16G RAM would be enough for a reasonably big dataset (\~ 40k cells, 10k genes) with multi-threading of 12 cpu cores. 
### Parameters 
(1) --maxloop  (optional but critical)  
Please set --maxloop 0 if you accept the clustering/trajectories from the prerun, it can dramatically cut-down the running time (from \~ hours -> \~ minutes).
The improvement varies a lot depending on specific application scenarios, in the above example, we got 10% improvement with spending 9x more running time. 
This parameter is optional, if you are not sure how to set it, just leave it as the default (or 0 as discussed above).       

(2) --logfc (optional) 
This specifies the log fold change cutoff for identifying and utilizing differential genes in the trajectories (and another condition:  student t-test p-value <0.05 is forced in the program)
Please set --logfc to a larger value (e.g., --logfc 1.5 -> 3x fold change) if you want to find and use "stringent" differential genes 
Please set --logfc to a smaller value (e.g., --logfc 0.6 --> 1.5x fold change) if you want to identify and use "possibly" more differential genes.    
We recommend a value better in a range 0.5 to 1.5.    This parameter is optional, if you are not sure how to set it, just leave it as the default.  

(3) --root (optional but critical)  
In a cell differentiation or disease progression process, you will always need to know what is the starting point (i.e., the root node) if we
assume the trajectory is tree-structured that is commonly accepted.   From the prerun results, we will find an umap clustering plot and PAGA trajectory as I mentioned
in the above example section.  Combining your prior knowledge (e.g, specific markers) with the prerun results, you will be able to make an informed choice of the root.  
You can check the expression for any specific gene (e.g, markers) with the provided "plotGene.py" script under the [utils](utils) folder. 

(4) --etfListFile (optional)  
By default, we provided a list of collected TFs in human and mouse, which is quite comprehensive with more than 1.2k TFs.  
If you are working in another species or you want to use your own list of TFs, please use this parameter. 
Please make sure the gene names are consistent with your input gene expression file.   
If you work in human and mouse, you can leave it as the default.  

(5) --llhCut (optional)  
the initial results will be iteratively refined by a PGM approach. In each loop, we will calculate a likelihood (of the cell assignment to the model). 
The likelihood will keep increasing and this parameter specifies when to stop.   By default, it's set as 0.05, which if the improvement is less than 5% of the original likelihood, it will stop.
Here, we recommend a value in a range from 0.01 (1%) to 0.1 (10%).  A smaller value would dramatically increase the running time while providing marginal improvement. 
A larger value would approximate --maxloop 0.  
  
# UTILS
We will keep provide custom scripts to better explore the model and predictions.  
All those utility scripts will be placed under [utils](utils) folder.   
All user contributed scripts are highly welcomed (best to provide an example test to prove its accuracy). 
You can add the script via the github or email me @ the address given below.  
__You credits__ will be __honored__ in the CREDITS section below.  

Current UTILS scripts:  
1)   [plotGene.py](utils/plotGene.py)  
This script plots the expression of input genes in different clusters   
2)   [buildInputExpressionFile.py](utils/buildInputExpressionFile.py)  
Build the tab-delimited expression.E file that required by the pre-run from 10x genomics output (mtx.tsv, barcodes.tsv, genes.tsv)

You use the -h command to look for the detailed usage of the scripts.  


# CREDITS
 
This software was developed by ZIV-system biology group @ Carnegie Mellon University.    
Developed by Jun Ding.  

Please cite our paper [Reconstructing differentiation networks and their regulation from time series single cell expression data](https://genome.cshlp.org/content/early/2018/01/09/gr.225979.117).  
scdiff2 uses many functions from scanpy (https://scanpy.readthedocs.io/en/stable/api/index.html) for the initialization.  Many thanks for the great tool.   
# LICENSE 
 
This software is under MIT license.  
see the LICENSE.txt file for details. 


# CONTACT 
jund  at cs.cmu.edu




                                 
                                 
                                 
                                 
                                 

                                                     
