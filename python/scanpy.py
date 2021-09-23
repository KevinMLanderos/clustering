#/usr/bin/env python

######################
# Clustering: Scanpy #
######################

import optparse
parser=optparse.OptionParser(usage="usage: %prog [options] --yaml config.yaml")
parser.add_option("-y", "--yaml", default='config.yaml',
                  help="Configuration file: Instructions in YAML format.",
                  metavar="FILE")
parser.add_option("-p", "--percent",
                  help="Percentage of variance.",
                  metavar="NUM")
parser.add_option("-n", "--n_comp", default=50,
                  help="Total number of components to explore.",
                  metavar="NUM")
parser.add_option("-c", "--chosen_comp",
                  help="Chosen components. Default: 20.",
                  metavar="NUM")
parser.add_option("-r", "--resolution",
                  help=("Resolution(s). Usually when running just the markers.  "
                        "It does not replace resolution in the YAML."),
                  metavar="NUM")
parser.add_option("-m", "--do_markers", action="store_true", default=True,
                  help="Find markers.",
                  metavar="True")
parser.add_option("-t", "--prefix",
                  help=("Prefix for this combination of parameters. Default is  "
                        "determined from the configuration file."),
                  metavar="True")
parser.add_option("-v", "--verbose", action="store_true", default=True,
                  help="Verbose: Show progress.",
                  metavar="True")

# Getting arguments from command line
(opt, args) = parser.parse_args()

import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import math
import os
import datetime

with open(opt.yaml, 'r') as file:
    try:
        config = yaml.safe_load(file)
    except yaml.YAMLError as exc:
        print(exc)

if opt.verbose:
  print("\n************ Vijay Lab - LJI")
  print("-------------------------------------")
  print("------------ Scanpy clustering")

#### Digesting parameters ###---------------------------------------------------
if config['dim_reduction']['base']['chosen_comp'] == None or opt.chosen_comp != None:
    config['dim_reduction']['base']['chosen_comp'] = opt.chosen_comp

if config['dim_reduction']['base']['chosen_comp'] == None:
    config['dim_reduction']['base']['chosen_comp'] = 0
#### ######### ########## ###---------------------------------------------------

#### Directories structure ####-------------------------------------------------
output_dir = config['output_dir'].rstrip("/") + "/" + config['project_name']
if not any(x in os.getcwd() for x in ["scratch", "beegfs"]):
    print("No scratch folder involved; careful about temp files...")
    if os.path.isdir(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)

print("Working in: %s" % os.getcwd())
if opt.prefix == None:
    global_prefix = "scanpy" # create_prefix
else:
    global_prefix = opt.prefix
#### ########### ######### ####-------------------------------------------------
### Log file ####---------------------------------------------------------------
print("Date and time:", datetime.datetime.now())
### ### #### ####---------------------------------------------------------------


# filters_file = ".seurat_filters.txt"
# these_filters = if(file.exists(filters_file)) readRDS(filters_file)
# the file that will store the analysis results
init_object_name = ".object_stem_" + global_prefix + ".h5ad"
these_filters == config['filtering']
if os.path.isfile(init_object_name):
    adata = sc.read_10x_mtx(
        config['input_expression'], # the directory with the `.mtx` file
        var_names='gene_symbols',   # use gene symbols for the variable names (variables-axis index)
        cache=True)                 # write a cache file for faster subsequent reading

    adata.var_names_make_unique()
    nSamples_min = math.ceil(config['filtering']['nSamples_expressed'] * adata.n_obs)
    sc.pp.filter_cells(adata, min_genes=0)
    sc.pp.filter_genes(adata, min_cells=nSamples_min)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        min_mean=config['variable_features']['mean.cutoff'][1],
        max_mean=config['variable_features']['mean.cutoff'][2],
        min_disp=config['variable_features']['dispersion.cutoff'][1])
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver='arpack')

    adata.write(init_object_name)
else:
    adata = sc.read(init_object_name)

tmpfile = open(".object_init_" + global_prefix, "w")
tmpfile.writelines(init_object_name)
file1.close() #to change file access modes

sc.pl.highly_variable_genes(adata).savefig('scanpy_hvg.png', dpi=300, bbox_inches='tight')

chosen_comp = config['dim_reduction']['base']['chosen_comp']
fname = global_prefix + "_elbow_" + chosen_comp + "of" + config['dim_reduction']['base']['n_comp'] + ".pdf"
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)

sc.tl.paga(adata)
sc.tl.umap(adata, init_pos='paga')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

adata.write(results_file)

adata = sc.read(results_file)
