"""
Aim: Snakemake workflow for clustering analyses
Date: Wednesday 11th of November 2020
Contact: Ciro Ramirez-Suastegui (ksuasteguic@gmail.com)
"""

import pandas as pd
import os
import sys

conf_file_loc = "config.yaml"

# Read config
configfile: conf_file_loc

# import yaml
# with open(conf_file_loc, 'r') as stream:
#     config = yaml.load(stream, Loader=yaml.FullLoader)

MEANS = config['variable_features']['mean.cutoff'][0]
PERCENTAGES = config['variable_features']['percent']
COMPONENTES = config['dim_reduction']['base']['chosen_comp']
RESOLUTIONS = config['resolution']
TOOLS = config['tool']
PIPELINE = config['pipeline'].rstrip("/")

try:
    df_exe = pd.read_csv(PIPELINE + "/data/tools.csv", index_col = 0)
except FileNotFoundError:
    sys.exit('Table of tools does not exists')

exec_report = PIPELINE + "/R/report.R"
exec_clustering = df_exe.loc[config['tool'], 'prefix_exec'] + " " + PIPELINE + "/" + df_exe.loc[config['tool'], 'script']

rule all:
    input:
        expand(".object_init_{tool}_mean{mean}_pct{percentage}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES),
        expand(".dr_{tool}_mean{mean}_pct{percentage}_pc{component}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
        expand(".markers_{tool}_mean{mean}_pct{percentage}_pc{component}_res{resolution}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES, resolution = RESOLUTIONS),
        expand(".{tool}_finished_components", tool = TOOLS, mean = MEANS),
        expand(".{tool}_finished_markers", tool = TOOLS, mean = MEANS)

rule init_object:
    input:
        conf_file_loc
    output:
        ".object_init_{tool}_mean{mean}_pct{percentage}"
    params:
        component = config['dim_reduction']['base']['n_comp']
    message: " --- Create initial object --- "
    shell:
        "{exec_clustering} -y {input} --percent {wildcards.percentage} --n_comp {params.component} "
        "--prefix {wildcards.tool}_mean{wildcards.mean}_pct{wildcards.percentage} --do_markers FALSE"

rule components:
    input:
        ".object_init_{tool}_mean{mean}_pct{percentage}"
    output:
        ".dr_{tool}_mean{mean}_pct{percentage}_pc{component}"
    message: " --- Branch number of components for clustering --- "
    shell:
        "{exec_clustering} -y {conf_file_loc} --chosen_comp {wildcards.component} "
        "--prefix {wildcards.tool}_mean{wildcards.mean}_pct{wildcards.percentage} --do_markers FALSE"

rule report_components:
    input:
        expand(".dr_{tool}_mean{mean}_pct{percentage}_pc{component}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES)
    output:
        ".{tool}_finished_components"
    message: " --- Creating report: components --- "
    shell:
        "Rscript {exec_report} --path ./ -m FALSE"

rule markers:
    input:
        ".dr_{tool}_mean{mean}_pct{percentage}_pc{component}"
    output:
        ".markers_{tool}_mean{mean}_pct{percentage}_pc{component}_res{resolution}"
    message: " --- Branch resolution for marker calculation --- "
    shell:
        "{exec_clustering} -y {conf_file_loc} --chosen_comp {wildcards.component} "
        "--prefix {wildcards.tool}_mean{wildcards.mean}_pct{wildcards.percentage} --resolution {wildcards.resolution}"

rule report_markers:
    input:
        expand(".markers_{tool}_mean{mean}_pct{percentage}_pc{component}_res{resolution}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES, resolution = RESOLUTIONS)
    output:
        ".{tool}_finished_markers"
    message: " --- Creating report: markers --- "
    shell:
        "Rscript {exec_report} --path ./ -c FALSE"
