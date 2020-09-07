# trajectory_inference
Both tutorial and basic pipeline for scRNA-seq trajectory analysis using dyno wrapper

## Info

* Article: https://www.nature.com/articles/s41587-019-0071-9
* Quick start: https://dynverse.org/users/2-quick_start/
* Methods: https://dynverse.org/reference/dynmethods/

![roadmap ti methods](https://dynverse.org/images/dyno/toolkit.png)

## Requirements

```{commandline}
r-argparse
dyno
dplyr
ggrepel
singularity>=3.0.1
snakemake
tidyverse
```

* dyno installation: https://dynverse.org/users/1-installation/
* singularity installation: https://anaconda.org/conda-forge/singularity

Before the first run, run this command in shell:

```{commandline}
export TAR=$(which tar)
```

## Guide

* `guides/guide.Rmd` -- markdown file which includes seurat object preparation using pbcm3k data and its
analysis using slingshot singularity container builded by dyno team
* `guides/trajectory analysis.html` -- the result of `guides/guide.Rmd` execution

## Pipeline

The basic pipeline for seurat object trajectory analysis was prepared using Snakemake

* tools -- tool(s) which you want to use
* assays -- target assay(s) -- RNA and (or) SCT
* sng_cache -- path to singularity cache directory

Pipeline execution:

```{commandline}
snakemake -j 2 --config tools=slingshot,tscan assays=SCT sng_cache=sin_cache
```