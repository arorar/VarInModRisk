rm(list=ls())

library(gridExtra)
library(ggplot2)
library(moments)

source("./Variance.R")
source("./Moments.R")
source("./Quantile.R")
source("./Risk.R")
source("./Plots.R")

# VAR
etl <- FALSE

# Gaussian large-sample
method <- "large-sample"; distribution <- "gaussian"
p1 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=etl)
# Gaussian Small Sample
method <- "small-sample"; distribution <- "gaussian"
p2 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=etl)

grid.arrange(p1, p2, nrow=2,ncol=1)

# t large-sample
method <- "large-sample"; distribution <- "t"
p3 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=etl)

# t Small Sample
method <- "small-sample"; distribution <- "t"
p4 <-plot.efficiency(distribution = distribution, method = method, 
                     plot = FALSE, etl=etl)

grid.arrange(p3, p4, nrow=2,ncol=1)

# ETL
etl <- TRUE

# Gaussian large-sample
method <- "large-sample"; distribution <- "gaussian"
p5 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=etl)
# Gaussian Small Sample
method <- "small-sample"; distribution <- "gaussian"
p6 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=etl)

grid.arrange(p5, p6, nrow=2,ncol=1)

# t large-sample
method <- "large-sample"; distribution <- "t"
p7 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=etl)

# t Small Sample
method <- "small-sample"; distribution <- "t"
p8 <-plot.efficiency(distribution = distribution, method = method, 
                     plot = FALSE, etl=etl)

grid.arrange(p7, p8, nrow=2,ncol=1)
