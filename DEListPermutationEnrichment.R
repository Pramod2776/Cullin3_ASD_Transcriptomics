if (!require(tidyverse, lib.loc = "tmp")) {
    install.packages("tidyverse", lib = "tmp", repos = "https://cran.rstudio.com")
}
if (!require(doParallel, lib.loc = "tmp")) {
    install.packages("doParallel", lib = "tmp", repos = "https://cran.rstudio.com")
}
library(tidyverse, lib.loc = "tmp")
library(doParallel, lib.loc = "tmp")

registerDoParallel(cores = detectCores() - 1)

permutation <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(permutation)) {
    permutation <- 1
}

deg_files <- list.files(path = "src", pattern = "*edgeR*", full.names = TRUE)
dep_files <- list.files(path = "src", pattern = "*limma*", full.names = TRUE)

three_period <- c("E17.5", "P7", "P35")
three_region <- c("CX", "CB", "HIP")
two_region <- c("CX", "CB")

region_deg_file <- lapply(setNames(nm = three_region), function(x){
  list.files(path = "./src",pattern = paste0(x,"_*edgeR*"), full.names = T)
})

period_deg_file <- lapply(setNames(nm = three_period), function(x){
  list.files(path = "./src",pattern = paste0("^",x,"_[A-Z]*_edgeR*"), full.names = T)
})

period_dep_file <- lapply(setNames(nm = three_period), function(x){
  list.files(path = "./src",pattern = paste0("^",x,"_[A-Z]*_CI"), full.names = T)
})

region_dep_file <- lapply(setNames(nm = two_region), function(x){
  list.files(path = "./src",pattern = paste0(x,"_CI"), full.names = T)
})
  

deg_regions <- mclapply(setNames(nm = names(region_deg_file)), function(f) {
    mclapply(region_deg_file[[f]], function(x){
      read_csv(x) %>%
        mutate(File = x) %>%
        filter(FDR <= 0.10)
    }) %>%
    bind_rows()
})

deg_periods <- mclapply(setNames(nm = names(period_deg_file)), function(f) {
  mclapply(period_deg_file[[f]], function(x){
    read_csv(x) %>%
      mutate(File = x) %>%
      filter(FDR <= 0.10)
  }) %>%
    bind_rows()
})

dep_regions <- mclapply(setNames(nm = names(region_dep_file)), function(f) {
    mclapply(region_dep_file[[f]], function(x){
      read_csv(x) %>%
        mutate(File = x) %>%
        filter(adj.P.Val <= 0.15)
    }) %>%
    bind_rows()
})

dep_periods <- mclapply(setNames(nm = names(period_dep_file)), function(f) {
  mclapply(period_dep_file[[f]], function(x){
    read_csv(x) %>%
      mutate(File = x) %>%
      filter(adj.P.Val <= 0.15)
  }) %>%
    bind_rows()
})

# gene_meta <- lapply(
#     list.files("src/gl_gc/", pattern = "*edgeR*", full.names = TRUE),
#     read_csv
# ) %>%
#     bind_rows() %>%
#     dplyr::select(Gene_name, percentage_gene_gc_content, gene_length) %>%
#     rename(gene_symbol = Gene_name) %>%
#     mutate(gene_symbol = toupper(gene_symbol)) %>%
#     bind_rows(
#         lapply(
#             list.files("src/gl_gc/", pattern = "*limma*", full.names = TRUE),
#             read_csv
#         ) %>%
#             bind_rows() %>%
#             rename(gene_symbol = GeneSymbol) %>%
#             mutate(gene_symbol = toupper(gene_symbol)) %>%
#             dplyr::select(gene_symbol, percentage_gene_gc_content, gene_length)
#     ) %>%
#     distinct()
# write_csv(gene_meta, "data/gene_length_and_gc.csv")
gene_metadata <- read_csv("data/gene_length_and_gc.csv") %>%
    filter(!is.na(gene_symbol))

#' Select a similar gene based on GC content and length
#' 
#' @param gene Input gene
#' @param gc Input gene GC content (percentage)
#' @param ln Input gene length (bp)
#' @param gc_min Lower-bound resampled GC content modifier (gc * gc_min)
#' @param gc_min Uppber-bound resampled GC content modifier (gc * gc_max)
#' @param ln_min Lower-bound resampled length modifier (ln * ln_min)
#' @param ln_min Upper-bound resampled length modifier (ln * ln_max)
#' @param all_genes Table of genes to resample from, must have cols
#'                  gene_symbol | percentage_gene_gc_content | gene_length
select_similar_genes <- function(
    gene, gc, ln, gc_min, gc_max, ln_min, ln_max, all_genes
) {
    minimum_gc <- gc - (gc * gc_min)
    maximum_gc <- gc + (gc * gc_max)
    minimum_ln <- ln - (ln * ln_min)
    maximum_ln <- ln + (ln * ln_max)
    resample_set <- all_genes %>%
        filter(percentage_gene_gc_content >= minimum_gc) %>%
        filter(percentage_gene_gc_content <= maximum_gc) %>%
        filter(gene_length >= minimum_ln) %>%
        filter(gene_length <= maximum_ln) %>%
        pull(gene_symbol)
    if (length(resample_set) == 0) {
        return(gene)
    } else {
        return(resample_set)
    }
}

# Do `RESAMPLING_ROUNDS` rounds of resampling
RESAMPLING_ROUNDS <- 10
deg_resample_region <- mclapply(
    setNames(nm = names(deg_regions)),
    function(x) {
        genes <- toupper(pull(deg_regions[[x]], Gene_name))
        lapply(1:RESAMPLING_ROUNDS, function(i) {
            message("Resampling, for ", x, ", iteration ", i)
            sapply(genes, function(g) { 
                gene_gc <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], percentage_gene_gc_content)[[1]]
                gene_ln <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], gene_length)[[1]]
                sample(select_similar_genes(g, gene_gc, gene_ln, 0.10, 0.10, 0.10, 0.10, gene_metadata), 1)
            })
        })
    }
)

deg_resample_period <- mclapply(
  setNames(nm = names(deg_periods)),
  function(x) {
    genes <- toupper(pull(deg_periods[[x]], Gene_name))
    lapply(1:RESAMPLING_ROUNDS, function(i) {
      message("Resampling, for ", x, ", iteration ", i)
      sapply(genes, function(g) { 
        gene_gc <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], percentage_gene_gc_content)[[1]]
        gene_ln <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], gene_length)[[1]]
        sample(select_similar_genes(g, gene_gc, gene_ln, 0.10, 0.10, 0.10, 0.10, gene_metadata), 1)
      })
    })
  }
)


dep_resample_region <- mclapply(
    setNames(nm = names(dep_regions)),
    function(x) {
        genes <- toupper(pull(dep_regions[[x]], GeneSymbol))
        lapply(1:RESAMPLING_ROUNDS, function(i) {
            message("Resampling, for ", x, ", iteration ", i)
            sapply(genes, function(g) { 
                gene_gc <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], percentage_gene_gc_content)[[1]]
                gene_ln <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], gene_length)[[1]]
                sample(select_similar_genes(g, gene_gc, gene_ln, 0.10, 0.10, 0.10, 0.10, gene_metadata), 1)
            })
        })
    }
)

dep_resample_period <- mclapply(
  setNames(nm = names(dep_periods)),
  function(x) {
    genes <- toupper(pull(dep_periods[[x]], GeneSymbol))
    lapply(1:RESAMPLING_ROUNDS, function(i) {
      message("Resampling, for ", x, ", iteration ", i)
      sapply(genes, function(g) { 
        gene_gc <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], percentage_gene_gc_content)[[1]]
        gene_ln <- pull(gene_metadata[gene_metadata$gene_symbol == g, ], gene_length)[[1]]
        sample(select_similar_genes(g, gene_gc, gene_ln, 0.10, 0.10, 0.10, 0.10, gene_metadata), 1)
      })
    })
  }
)


# Calculate overlaps per list, per round of resampling
gene_lists <- lapply(
    setNames(nm = readxl::excel_sheets('src/ASDRelevantGeneListsFromLiterature.xlsx')[-13]),
    function(sheet) {
        readxl::read_xlsx(
            'src/ASDRelevantGeneListsFromLiterature.xlsx',
            sheet = sheet
        )[, 1] %>%
            distinct() %>%
            pull(gene_symbol)
    }
)
gene_lists <- gene_lists[!names(gene_lists) %in% c(
    "EichlerDNM_LGD_MIS_ASD_BD_SCZ", 
    "EichlerDNM_LGD_ASD_BD_SCZ",
    "VulnerableASD",
    "ASDSanders65"
)]

null_distribution_deg_region <- foreach (i=1:length(names(deg_regions))) %dopar% {
    deg_fn <- names(deg_regions)[[i]]
    this_file_permutation <- tibble()
    for (j in 1:length(deg_resample_region[[deg_fn]])) {
        this_permutation <- tibble()
        for (gl in names(gene_lists)) {
            overlap <- length(intersect(deg_resample_region[[deg_fn]][[j]], gene_lists[[gl]]))
            this_permutation <- bind_rows(
                this_permutation,
                tibble(FILE = deg_fn, RESAMPLE = j, GENELIST = gl, OVERLAP = overlap)
            )
        }
        this_file_permutation <- bind_rows(this_file_permutation, this_permutation)
    }
    return(this_file_permutation)
} %>%
    bind_rows()
write_csv(
    null_distribution_deg_region,
    paste0("data/DE_NULL_DISTRIBUTIONS/null_distribution_deg_PERMUTATION_", permutation, "_RegionWise.csv")
)

null_distribution_deg_period <- foreach (i=1:length(names(deg_periods))) %dopar% {
  deg_fn <- names(deg_periods)[[i]]
  this_file_permutation <- tibble()
  for (j in 1:length(deg_resample_period[[deg_fn]])) {
    this_permutation <- tibble()
    for (gl in names(gene_lists)) {
      overlap <- length(intersect(deg_resample_period[[deg_fn]][[j]], gene_lists[[gl]]))
      this_permutation <- bind_rows(
        this_permutation,
        tibble(FILE = deg_fn, RESAMPLE = j, GENELIST = gl, OVERLAP = overlap)
      )
    }
    this_file_permutation <- bind_rows(this_file_permutation, this_permutation)
  }
  return(this_file_permutation)
} %>%
  bind_rows()
write_csv(
  null_distribution_deg_period,
  paste0("data/DE_NULL_DISTRIBUTIONS/null_distribution_deg_PERMUTATION_", permutation, "_PeriodWise.csv")
)

null_distribution_dep_region <- foreach (i=1:length(names(dep_regions))) %dopar% {
    dep_fn <- names(dep_regions)[[i]]
    this_file_permutation <- tibble()
    for (j in 1:length(dep_resample_region[[dep_fn]])) {
        this_permutation <- tibble()
        for (gl in names(gene_lists)) {
            overlap <- length(intersect(dep_resample_region[[dep_fn]][[j]], gene_lists[[gl]]))
            this_permutation <- bind_rows(
                this_permutation,
                tibble(FILE = dep_fn, RESAMPLE = j, GENELIST = gl, OVERLAP = overlap)
            )
        }
        this_file_permutation <- bind_rows(this_file_permutation, this_permutation)
    }
    return(this_file_permutation)
} %>%
    bind_rows()
write_csv(
    null_distribution_dep,
    paste0("data/DE_NULL_DISTRIBUTIONS/null_distribution_dep_PERMUTATION_", permutation, "_RegionWise.csv")
)


null_distribution_dep_period <- foreach (i=1:length(names(dep_periods))) %dopar% {
  dep_fn <- names(dep_periods)[[i]]
  this_file_permutation <- tibble()
  for (j in 1:length(dep_resample_period[[dep_fn]])) {
    this_permutation <- tibble()
    for (gl in names(gene_lists)) {
      overlap <- length(intersect(dep_resample_period[[dep_fn]][[j]], gene_lists[[gl]]))
      this_permutation <- bind_rows(
        this_permutation,
        tibble(FILE = dep_fn, RESAMPLE = j, GENELIST = gl, OVERLAP = overlap)
      )
    }
    this_file_permutation <- bind_rows(this_file_permutation, this_permutation)
  }
  return(this_file_permutation)
} %>%
  bind_rows()
write_csv(
  null_distribution_dep,
  paste0("data/DE_NULL_DISTRIBUTIONS/null_distribution_dep_PERMUTATION_", permutation, "_PeriodWise.csv")
)

message("Null distributions calculations complete")