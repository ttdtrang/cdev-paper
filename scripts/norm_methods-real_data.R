devtools::load_all('./')
library(magrittr)
library(parallel)


NCORES = 6
load('data/data_container.Rdata')

experimental_conditions = c(
    'rbm1' = 'Organ',
    'rbm2' = 'Organ',
    'rnor1' = 'Chemical',
    'rnor2' = 'Chemical',
    'mpreg1a' = 'BrainRegion',
    'mpreg1b' = 'Stage',
    'mpreg2a' = 'BrainRegion',
    'mpreg2b' = 'Stage',
    'dmelAC' = 'Environment',
    'dmelAV' = 'Environment',
    'sarcoma' = 'Tissue',
    'lymphoma' = 'Treatment',
    'yfv' = 'Infection',
    'xtro1m' = 'Time.hpf.',
    'xtro2' = 'Time.hpf.'
)

start = proc.time()
env.all_methods = new.env()

# Definitions of various normalization methods

norm_methods = list('raw' = raw,
                    'TC' = tc,
                    'UQ' = uq.v2,
                    'TMM' = tmm.v2,
                    'DESeq' = deseq.v1,
                    'DEGES/TMM' = deges.3,
                    'PoissonSeq' = pseq.v1,
                    'Ground-truth' = NA)
norm_methods$raw <- function(X, group=NA) {
    return(X)
}

# Run

for (data_prefix in names(experimental_conditions)) {

    condition = experimental_conditions[data_prefix]
    X = get(paste0(data_prefix,'.cnt'), envir = data.container)
    group = get(paste0(data_prefix,'.pheno'), envir = data.container)[[condition]]
    REFS.idx = get(paste0(data_prefix,'.refs.idx'), envir = data.container)

    norm_methods[['Ground-truth']] <- function(X, group = NA) {
        return(normalize.by.refs(X, REFS.idx))
    }

    # norm_methods$loess <- function(X, group = NA) {
    #     return(affy.loess(X, REFS.idx))
    # }
    #
    # norm_methods$RUVg <- function(X, group = NA) {
    #     return(RUVSeq::RUVg(log2(X+1),cIdx = REFS.idx, k = 1, isLog = TRUE))
    # }
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

    assign(paste0(data_prefix, '.normResults'), x.normResults, envir = env.all_methods)
    assign(paste0(data_prefix, '.cdevList'), x.cdevList, envir = env.all_methods)
    assign(paste0(data_prefix, '.cdevMatrix'), x.cdevMatrix, envir = env.all_methods)

}

dt = proc.time() - start
print(sprintf("Time to compute: %s", dt['elapsed']))
saveRDS(env.all_methods, file = 'real_data_normalized.RDS')

