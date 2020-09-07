suppressMessages(library(argparse))
suppressMessages(library(babelwhale))
suppressMessages(library(dyno))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))


parser <-
  ArgumentParser(description = 'Trajectory analysis using dyno wrapper')
parser$add_argument('--assay',
                    type = "character",
                    help = 'Assay name: RNA or SCT')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to the Seurat object in RData format')
parser$add_argument('--ti_tool',
                    type = "character",
                    help = 'Name of the trajectory inference tool available at dyno wrapper')
parser$add_argument('--sng_cache',
                    type = "character",
                    help = 'Path to the cahce dir for singularity container')
## SET VARIABLES

args <- parser$parse_args()

assay <- args$assay
ti_tool <- args$ti_tool

print(args)

## CONFIG SINGULARITY

babelwhale::set_default_config(babelwhale::create_singularity_config(cache_dir = args$sng_cache))

## DEFINE SOME FUCNTIONS

get_data <- function(path, assay) {
  obj <- get(load(path))
  wrap_data <- wrap_expression(counts = t(as.matrix(obj@assays[[assay]]@counts)),
                               expression = t(as.matrix(obj@assays[[assay]]@data)))
  wrap_data <- add_grouping(wrap_data,
                            obj$seurat_clusters)
  wrap_data
}

get_branching_point <- function(model, obj) {
  branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% dplyr::first()
  branch_point <- calculate_branching_point_feature_importance(model,
                                                               expression_source=obj$expression,
                                                               milestones_oi = branching_milestone)
  branch_point
}

get_hmap <- function(to, branch_feature_importance, model, wrap_data, topn=30) {
  features <- branch_feature_importance %>%
    filter(to == to) %>%
    top_n(topn, importance) %>%
    pull(feature_id)
  plot_heatmap(model, expression_source = wrap_data, features_oi = unique(features))+
    ggtitle(sprintf('pbmc, to: %s', to))
  ggsave(sprintf('%s/%s/plots/to-%s.png', ti_tool, assay, to))
}

save_tables <- function(model, ti_tool, assay, branch_feature_importance, branch_point) {
  write.table(model$progressions,
              file=sprintf('%s/%s/tables/progression.csv', ti_tool, assay), quote = F, row.names = F, sep=',')
  write.table(branch_feature_importance,
              file=sprintf('%s/%s/tables/branch_feature_importance.csv', ti_tool, assay), quote = F, row.names = F, sep=',')
  write.table(branch_point,
              file=sprintf('%s/%s/tables/branch_point.csv', ti_tool, assay), quote = F, row.names = F, sep=',')
}

## PREPARE THE OBJECT

obj <- get_data(args$data, assay)

## INFER TRAJECTORY

set.seed(1)

model <- infer_trajectory(obj, method = args$ti_tool, verbose = TRUE)

## UMAP

dimred <- dyndimred::dimred_umap(obj$expression)

## CALCULATE SOME FEATURES: FEATURE IMPORTANCE, BRANCHING POINT

branch_feature_importance <- calculate_branch_feature_importance(model, expression_source=obj$expression)
branch_point <- get_branching_point(model, obj)

## VISUALIZATION

plot_dimred(
  model,
  dimred = dimred,
  expression_source = obj$expression,
  color_density = "grouping",
  alpha_cells = 0.7,
  grouping = obj$grouping,
  label_milestones = FALSE,
)+
  theme(aspect.ratio = 1, legend.position = "none")
ggsave(sprintf('%s/%s/plots/trajectory.png', ti_tool, assay))


sapply(unique(branch_feature_importance$to), function(x) get_hmap(x, branch_feature_importance, model, obj))

## SAVE OUTPUT: TABLES

save_tables(model, ti_tool, assay, branch_feature_importance, branch_point)

## SAVE OUTPUT: RDATA

save(list = c('args', 'obj', 'model', 'dimred'),
     file = sprintf('%s/%s/rdata/pbmc.RData', ti_tool, assay))


