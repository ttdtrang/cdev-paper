devtools::load_all('./')
library(magrittr)
library(parallel)

NSAMPLING_FULL = 30
NSAMPLING_SUBSET = 30
NCORES = 6
load('data_container.Rdata')
# DATASET_PREFIXES = c('rat1', 'rat2', 'rbm1', 'rbm2', 'rnor1', 'rnor2',
#                      'mpreg1a', 'mpreg1b', 'mpreg2a', 'mpreg2b',
#                      'dmelA', 'dmelAC', 'dmelAV',
#                      'sarcoma', 'lymphoma', 'yfv',
#                      'xtro1', 'xtro1m', 'xtro2', 'xtroAp', 'xtroAr', 'xtroB')
DATASET_PREFIXES = c('dmelA', 'dmelAC', 'dmelAV')

start = proc.time()
env.randrefs_norm = new.env()
for (data_prefix in DATASET_PREFIXES) {
    nrefs = length(get(paste0(data_prefix, '.refs.idx'), envir = data.container))
    x.cnt = get(paste0(data_prefix, '.cnt'), envir = data.container)

    is.nonZero = (apply(x.cnt,1, min) > 0)
    is.ERCC = grepl('^ERCC', rownames(x.cnt))
    lapply(1:NSAMPLING_FULL, function(j) {
        rand.refs = sample(which(is.nonZero & (!is.ERCC)), size = nrefs)
        x.randref = normalize.by.refs(x.cnt, ref.idx = rand.refs)

        fullset_avgCV = mean(apply(x.cnt[rand.refs,,drop=FALSE], 1, cv))
        fullset_totalVar = sum(apply(x.cnt[rand.refs,,drop=FALSE], 1, var))
        # n_combinations = min(10, choose(nrefs, floor(nrefs/2)))
        # tryCatch({
            mclapply(1:NSAMPLING_SUBSET, mc.cores = NCORES, function(i) {
                rand_subset = sample(rand.refs, size = floor(nrefs/2))
                x.rand_subset = normalize.by.refs(x.cnt, ref.idx = rand_subset)
                data.frame(
                    'random_fullset_index' = j,
                    'random_subset_index' = i,
                    'cdev' = cdev(x.rand_subset, x.randref),
                    'fullset_avgCV' = fullset_avgCV,
                    'fullset_totalVar' = fullset_totalVar,
                    'subset_CV' = cv(x.cnt[rand_subset,,drop=FALSE]),
                    'subset_avgCV' = mean(apply(x.cnt[rand_subset,,drop=FALSE], 1, cv)),
                    'subset_totalVar' = sum(apply(x.cnt[rand_subset,,drop=FALSE], 1, var)),
                    'subset_avgVar' = mean(apply(x.cnt[rand_subset,,drop=FALSE], 1, var))
                )
            }) %>%
                do.call(rbind, .)
        # }, error = function(e) {
            # message(e)
        # })
    }) %>%
        do.call(rbind, .) %>%
        assign(data_prefix, value = ., envir = env.randrefs_norm)
}

saveRDS(env.randrefs_norm, file = 'cdev_random_vs_random.2.RDS')

dt = proc.time() - start
print(sprintf("Time to compute: %s", dt['elapsed']))
