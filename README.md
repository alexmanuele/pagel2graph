# pagel2graph
Tools for converting Pagel results to GraphML format and subsequently filter the GraphML.

## Requirements
Requires pandas and networkx.

## Installation
Download the repository
```
git clone https://github.com/alexmanuele/pagel2graph.git
```
Suggested: Create a conda environment to manage dependancies
```
conda create -n pagel2graph networkx
conda activate pagel2graph
conda install pandas
```

## Usage
Expects as input a GraphML file where all nodes contain an "lr" attribute referring to likelihood ratio and a 'p' attribute referring to statistical p-value, both having been calculated from pagel.
Given a GraphML file and a node of interest, filter the GraphML file to contain the node and any neighbors satisfying edge filtering criteria up to a specified depth.

```
python filter_graphml.py \\
 -i input_file_path \\
 -n name of node of interest  \\
 -d degee of graph traversal (int) \\
 -lr minimum likelihood ratio. Values lower than this will be filtered \\
 -p maximum p value. Values higher than this will be filtered \\
 -o Output file name. \\
 ```
