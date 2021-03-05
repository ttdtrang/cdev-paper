devtools::load_all('./')
library(magrittr)
library(parallel)

NCORES = 6
NSAMPLING = 60
load('data_container.Rdata')
# DATASET_PREFIXES = c('rat1', 'rat2', 'rbm1', 'rbm2', 'rnor1', 'rnor2',
#                      'mpreg1a', 'mpreg1b', 'mpreg2a', 'mpreg2b',
#                      'dmelA', 'dmelAC', 'dmelAV',
#                      'sarcoma', 'lymphoma', 'yfv',
#                      'xtro1', 'xtro1m', 'xtro2', 'xtroAp', 'xtroAr', 'xtroB')
DATASET_PREFIXES = c('dmelA', 'dmelAC', 'dmelAV')


start = proc.time()
env.rand_norm = new.env()
for (data_prefix in DATASET_PREFIXES) {
    nrefs = length(get(paste0(data_prefix, '.refs.idx'), envir = data.container))
    x.normed = get(paste0(data_prefix, '.normed'), envir = data.container)
    x.cnt = get(paste0(data_prefix, '.cnt'), envir = data.container)
    is.nonZero = (apply(x.cnt,1, min) > 0)
    tryCatch({
        mclapply(1:NSAMPLING, mc.cores = NCORES, FUN = function(i) {
            rand.refs = sample(which(is.nonZero), size = nrefs)
            x.rand = normalize.by.refs(x.cnt, ref.idx = rand.refs)
            refs_cv_stats = fivenum(apply(x.cnt[rand.refs,,drop=FALSE], 1, cv))
            data.frame(
                'sampling_index' = i,
                'cdev' = cdev(x.rand, x.normed),
                'n_refs' = nrefs,
                'refs_totalVar' = sum(apply(x.cnt[rand.refs,,drop=FALSE], 1, var)),
                'refs_avgVar' = mean(apply(x.cnt[rand.refs,,drop=FALSE], 1, var)),
                'refs_avgCV' = mean(apply(x.cnt[rand.refs,,drop=FALSE], 1, var)),
                'refs_cv_min' = refs_cv_stats[1],
                'refs_cv_lq' = refs_cv_stats[2],
                'refs_cv_median' = refs_cv_stats[3],
                'refs_cv_uq' = refs_cv_stats[4],
                'refs_cv_max' = refs_cv_stats[5]
            )
        }) %>%
            do.call(rbind, .) %>%
            assign(data_prefix, value = ., envir = env.rand_norm)
    }, error = function(e) {
        message(e)
    })

}

dt = proc.time() - start
print(sprintf("Time to compute: %s", dt['elapsed']))
saveRDS(env.rand_norm, file = 'cdev_random_vs_ground.2.RDS')
