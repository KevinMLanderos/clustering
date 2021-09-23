## Clustering

All the information goes into the configuration file (YAML format). There is an example (config.yaml) with comments regarding the files' format.

There are two main R scripts:
- _demultiplexing_cells.R_, sets and runs the Cell Ranger routines count and vdj for each sample in 'fastqs_dir' (and as indicated in 'samples').
- _aggregate.R_, according to 'aggregation'.

The files you need to prepare are or how to prepare:
1. Sample metadata ('metadata' in the YAML file).
	- This is optional but it is highly recommended that you have this information at hand so the scripts can produce a full more complete report of the (for now categorical) variables of interest. It accepts more than one file.
2. Filtering has two options. First, use 'subset: {expr[ession]: "filter expression"}' to indicate a logical expression. Second, create a new item for each column you want to filter by, e. g., 'subset: {doublet: "no", expr: "logical expression"}'. Also, you can combine both options.

NOTES:
1. At the moment only Seurat is included as an option.

### Install
Clone this repository (your ~/bin folder is a good place).
```
git clone https://git .. /clustering.git [copy clone URL at the top]
cd clustering
```

### Run the pipeline
After you've added the necessary information to the YAML file you can call the pipeline.
```
run.sh -y /path/to/project/config.yaml
```

You might benefit from having an environment for the pipeline. And you can indicate the name in the YAML
file `environment: custom_name` at any point in the file. In this case it's called clustering
First, you need minconda:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh 
```

For changes to take effect, close and re-open your current shell.
Now, you can create the invironment:

```
conda create --name clustering python=3.8
conda activate clustering
pip install -r requirements.txt
conda deactivate
```

Finally, you install snakemake: `conda install snakemake`
