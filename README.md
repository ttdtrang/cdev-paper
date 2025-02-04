# cdev: A ground-truth based measure to evaluate RNA-seq normalization performance

R codes accompanying paper _cdev: A ground-truth based measure to evaluate RNA-seq normalization performance_

R notebooks:
 
* [examine-real-datasets.Rmd](./notes/examine-real-datasets.Rmd): To examine the behavior of spike-ins in various RNA-seq data sets, and collect usable data sets into a single file 
* [cdev-upstream-downstream.Rmd](./notes/cdev-upstream-downstream.Rmd): cdev in relation to upstream processing quality and downstream performance, i.e. DE analysis
* [cdev-real-data.Rmd](./notes/cdev-real-data.Rmd): behavior of cdev in real data sets
* [cdev-application.Rmd](./notes/cdev-application.Rmd): application of cdev to compare common normalization methods on simulation and real data

R scripts:
* [cdev-ground_vs_ground.R](./scripts/cdev-ground_vs_ground.R)
* [cdev-random_vs_ground.R](./scripts/cdev-random_vs_ground.R)
* [cdev-random_vs_random.R](./scripts/cdev-random_vs_random.R)
* [norm_methods-real_data.R](./scripts/norm_methods-real_data.R)
* [norm_methods-simulation.R](./scripts/norm_methods-simulation.R)
* 
