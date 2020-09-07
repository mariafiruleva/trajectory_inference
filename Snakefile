assays = config['assays'].split(',')
ti_tools = config['tools'].split(',')

rule all:
    input: expand("{ti_tool}/{assay}/rdata/pbmc.RData", ti_tool=ti_tools, assay=assays)

rule ti_analysis:
    input: data="data/pbmc.RData"
    output: rda="{ti_tool}/{assay}/rdata/pbmc.RData",
          plt="{ti_tool}/{assay}/plots/trajectory.png",
          tbls="{ti_tool}/{assay}/tables/progression.csv"
    benchmark: "benchmarks/{ti_tool}_{assay}.txt"
    log: "logs/{ti_tool}_{assay}.log"
    params: sng_cache=config['sng_cache']
    shell:
         """
         /scratch/opt/R/3.6.0/bin/Rscript scripts/ti_analysis.R --data {input.data} --assay {wildcards.assay} --ti_tool {wildcards.ti_tool} \
         --sng_cache {params.sng_cache} 2> {log}
         """
