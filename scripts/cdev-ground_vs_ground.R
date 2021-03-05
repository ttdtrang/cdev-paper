devtools::load_all('./')
library(magrittr)
library(parallel)


NCORES = 6
NSAMPLING = 30
load('data_container.Rdata')
# DATASET_PREFIXES = c('rat1', 'rat2', 'rbm1', 'rbm2', 'rnor1', 'rnor2',
#                      'mpreg1a', 'mpreg1b', 'mpreg2a', 'mpreg2b',
#                      'dmelA', 'dmelAC', 'dmelAV',
#                      'sarcoma', 'lymphoma', 'yfv',
#                      'xtro1', 'xtro1m', 'xtro2', 'xtroAp', 'xtroAr', 'xtroB')
DATASET_PREFIXES = c('dmelA', 'dmelAC', 'dmelAV')

start = proc.time()
env.subset_ground = new.env()
for (data_prefix in DATASET_PREFIXES) {
    true_refs = get(paste0(data_prefix, '.refs.idx'), envir = data.container)
    nrefs = length(true_refs)
    x.normed = get(paste0(data_prefix, '.normed'), envir = data.container)
    x.cnt = get(paste0(data_prefix, '.cnt'), envir = data.container)

    tryCatch({
        lapply(1:floor(nrefs/2), FUN = function(i) {
            mclapply(1:min(NSAMPLING,choose(nrefs,i)), mc.cores = NCORES, FUN = function(j) {
                rand.refs = sample(true_refs, size = i)
                x.rand = normalize.by.refs(x.cnt, ref.idx = rand.refs)
                refs_cv_stats = fivenum(apply(x.cnt[rand.refs,,drop=FALSE], 1, cv))
                data.frame(
                    'sampling_index' = j,
                    'cdev' = cdev(x.rand, x.normed),
                    'n_refs' = i,
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
                do.call(rbind, .)
        }) %>%
            do.call(rbind, .) %>%
            assign(data_prefix, value = ., envir = env.subset_ground)
    }, error = function(e) {
        message(e)
    })

}

dt = proc.time() - start
print(sprintf("Time to compute: %s", dt['elapsed']))
saveRDS(env.subset_ground, file = 'cdev_ground_vs_ground.2.RDS')
