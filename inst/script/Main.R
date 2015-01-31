rm(list=ls())

# Gaussian large-sample
method <- "large-sample"; distribution <- "gaussian"
p1 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=FALSE)
p2 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=TRUE)
grid.arrange(p1, p2, nrow=2,ncol=1)

# t large-sample
method <- "large-sample"; distribution <- "t"
p3 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=FALSE)
p4 <-plot.efficiency(distribution = distribution, method = method, 
                     plot = FALSE, etl=TRUE)

grid.arrange(p3, p4, nrow=2,ncol=1)

# Gaussian small-sample
method <- "small-sample"; distribution <- "gaussian"
p5 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=FALSE)
p6 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=TRUE)

grid.arrange(p5, p6, nrow=2,ncol=1)

# t small-sample
method <- "small-sample"; distribution <- "t"
p7 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=FALSE)

p8 <-plot.efficiency(distribution = distribution, method = method, 
                     plot = FALSE, etl=TRUE)

grid.arrange(p7, p8, nrow=2,ncol=1)
