# pagel2graph
Tools for converting Pagel results to GraphML format and subsequently filter the GraphML.

## Installation
Installation not currently supported for Windows.

Download the repository
```
git clone https://github.com/alexmanuele/pagel2graph.git
```

Suggested: Create a conda environment to manage dependancies. This requires Anaconda or Miniconda.

### Linux
```
conda env create -f environment.yml
```
### MacOS
Conda environment files created on Ubuntu machine don't work properly on Mac.
For now, to install:
```
conda create -n pagel2graph networkx
conda activate pagel2graph
conda install -c conda-forge dash
conda install -c conda-forge dash-bootstrap-components
conda install -c conda-forge dash_cytoscape
```


## Usage
This repository contains two utilities; a graphML filtering script and an interactive Dash application.

### Filtering
Expects as input a GraphML file where all nodes contain an "lr" attribute referring to likelihood ratio and a 'p' attribute referring to statistical p-value, both having been calculated from pagel.
Given a GraphML file and a node of interest, filter the GraphML file to contain the node and any neighbors satisfying edge filtering criteria up to a specified depth.

```
conda activate pagel2graph

python filter_graphml.py \\
 -i input_file_path \\
 -n name of node of interest  \\
 -d degee of graph traversal (int) \\
 -lr minimum likelihood ratio. Values lower than this will be filtered \\
 -p maximum p value. Values higher than this will be filtered \\
 -o Output file name. \\
 ```

### Dash Application

##### Configuration
The current iteration of this software uses hardcoded file paths. Please follow the below instructions to configure the data properly:

```
# Make sure you're in the right directory
cd pagel2graph
# Make a 'data' directory
mkdir data
```
The data directory must contain the following files:
```
efaecium_profile_LR_rerunNA.csv
efaecium_profile_pval_rerunNA.csv
pagel_LR_featureVsHabitat.csv
pagel_pvalue_featureVsHabitat.csv
pagel_results_as_network_updated.graphml
```
##### Usage
Once your data is configured properly, usage is simple.
Open a terminal and activate your conda environment:

`
conda activate pagel2graph
`
Then, type this command into your terminal from the root directory of the codebase:

`python app.py`

This will launch the Dash server. Next, simply open your favorite web browser and navigate to http://localhost:8050/ .
