"""
Aim: Snakemake workflow for clustering analyses
Date: Wednesday 11th of November 2020
Contact: Ciro Ramirez-Suastegui (ksuasteguic@gmail.com)
"""

import pandas as pd
import os, sys, re

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
EXECS = config['exec']
try:
    EXEC_R = list(filter(re.compile("Rscript").match, [EXECS] if type(EXECS)==str else EXECS))[0]
except IndexError:
    EXEC_R = "Rscript"
SCRIPTS = config['script']
EXEC_REPORT = config['pipeline'].rstrip("/") + "/R/report.R"

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
        "{EXECS} {SCRIPTS} -y {input} --percent {wildcards.percentage} --n_comp {params.component} "
        "--prefix {wildcards.tool}_mean{wildcards.mean}_pct{wildcards.percentage} --do_markers FALSE"

rule components:
    input:
        ".object_init_{tool}_mean{mean}_pct{percentage}"
    output:
        ".dr_{tool}_mean{mean}_pct{percentage}_pc{component}"
    message: " --- Branch number of components for clustering --- "
    shell:
        "{EXECS} {SCRIPTS} -y {conf_file_loc} --chosen_comp {wildcards.component} "
        "--prefix {wildcards.tool}_mean{wildcards.mean}_pct{wildcards.percentage} --do_markers FALSE"

rule report_components:
    input:
        expand(".dr_{tool}_mean{mean}_pct{percentage}_pc{component}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES)
    output:
        ".{tool}_finished_components"
    message: " --- Creating report: components --- "
    shell:
        "{EXEC_R} {EXEC_REPORT} --path ./ -m FALSE"

rule markers:
    input:
        ".dr_{tool}_mean{mean}_pct{percentage}_pc{component}"
    output:
        ".markers_{tool}_mean{mean}_pct{percentage}_pc{component}_res{resolution}"
    message: " --- Branch resolution for marker calculation --- "
    shell:
        "{EXECS} {SCRIPTS} -y {conf_file_loc} --chosen_comp {wildcards.component} "
        "--prefix {wildcards.tool}_mean{wildcards.mean}_pct{wildcards.percentage} --resolution {wildcards.resolution}"

rule report_markers:
    input:
        expand(".markers_{tool}_mean{mean}_pct{percentage}_pc{component}_res{resolution}", tool = TOOLS, mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES, resolution = RESOLUTIONS)
    output:
        ".{tool}_finished_markers"
    message: " --- Creating report: markers --- "
    shell:
        "{EXEC_R} {EXEC_REPORT} --path ./ -c FALSE"
