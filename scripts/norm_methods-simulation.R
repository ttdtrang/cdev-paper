devtools::load_all('./')
library(magrittr)
library(parallel)


NCORES = 1

## Varying DE percentage
NSAMPLES = 40
nGene = 10000
percentDEs = c(1:9/10)

sims_rp <- lapply(percentDEs, function(percentDE) {
    simulate.counts.varyingfc(ngenes = nGene, nlibs1 = NSAMPLES %/% 2, nlibs2= NSAMPLES - (NSAMPLES %/% 2), percentDE = percentDE)
}) %>%
    set_names(as.character(percentDEs))

sims_rp.ref = lapply(percentDEs, function(percentDE) {
    sim = sims_rp[[as.character(percentDE)]]
    normalize.by.refs(sim$counts, which(sim$differential == 0))
}) %>%
    set_names(as.character(percentDEs))

start = proc.time()

# Definitions of various normalization methods

norm_methods = list('raw' = raw,
                    'TC' = tc,
                    'UQ' = uq.v2,
                    'TMM' = tmm.v2,
                    'DESeq' = deseq.v1,
                    'DEGES/TMM' = deges.3,
                    'PoissonSeq' = pseq.v1,
                    # 'loess' = NA,
                    # 'RUVg' = NA,
                    'Oracle' = NA)
norm_methods$raw <- function(X, group=NA) {
    return(X)
}

# Run

dat.all_methods = lapply(percentDEs, function(percentDE) {
    sim = sims_rp[[as.character(percentDE)]]
    X = sim$counts
    group = sim$conditions$condition
    REFS.idx = which(sim$differential == 0)

    norm_methods$Oracle <- function(X, group = NA) {
        return(normalize.by.refs(X, REFS.idx))
    }

    x.normResults =
        mclapply(1:length(norm_methods),
                 mc.cores = NCORES,
                 # mc.cores = min(length(norm_methods), detectCores() - 1),
                 FUN = function(i) {
                     do.call(norm_methods[[i]], list(X, group = group))
                 }) %>%
        set_names(names(norm_methods))

    coords = combn(1:length(norm_methods), 2)
    x.cdevList = mclapply(1:ncol(coords),
                          mc.cores = NCORES,
                          # FUN.VALUE = 0.,
                          FUN = function (k) {
                              i = names(norm_methods[coords[1,k]])
                              j = names(norm_methods[coords[2,k]])
                              cdev(x.normResults[[i]], x.normResults[[j]])
                          }) %>%
        unlist()
    x.cdevMatrix = Matrix::sparseMatrix(i = coords[1,], j = coords[2,], x = x.cdevList,
                                        dimnames = list(names(norm_methods), names(norm_methods)),
                                        symmetric = TRUE)
    diag(x.cdevMatrix) = 1

    x.result = list('normResults' = x.normResults,
                    'cdevMatrix' = x.cdevMatrix,
                    'refs.idx' = REFS.idx)

}) %>%
    set_names(as.character(percentDEs))

dt = proc.time() - start
print(sprintf("Time to compute: %s", dt['elapsed']))
saveRDS(dat.all_methods, file = 'simulations_normalized.RDS')

