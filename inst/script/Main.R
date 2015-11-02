rm(list=ls())

# # Gaussian large-sample
method <- "large-sample";distribution <- "gaussian" 

# Gaussian large-sample VaR
p1 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=FALSE)

# Gaussian large-sample ETL
p2 <- plot.efficiency(distribution = distribution, method = method ,
                      plot = FALSE, etl=TRUE)

grid.arrange(p1, p2, nrow=2,ncol=1)

# t large-sample
method <- "large-sample";distribution <- "t" 

# student-t large-sample VaR
p1 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=FALSE)

# student-t large-sample ETL
p2 <- plot.efficiency(distribution = distribution, method = method, 
                      plot = FALSE, etl=TRUE)

grid.arrange(p1, p2, nrow=2,ncol=1)

# Gaussian small-sample
method <- "small-sample";distribution <- "gaussian" 

# Gaussian small-sample VaR
p1 <- plot.efficiency(distribution = distribution, method = method, 
                      size = c(100,250,500), plot = FALSE, etl=FALSE, 
                      replicates = 10000)

# Gaussian small-sample ETL
p2 <- plot.efficiency(distribution = distribution, 
                      method = method, size = c(100,250,500), plot = FALSE, 
                      etl=TRUE, replicates = 10000)

grid.arrange(p1, p2, nrow=2,ncol=1)

# t small-sample
method <- "small-sample";distribution <- "t" 

# Gaussian small-sample VaR
p1 <- plot.efficiency(distribution = distribution, method = method, 
                      size = c(100,250,500), plot = FALSE, etl=FALSE, 
                      replicates = 10000)

# Gaussian small-sample ETL
p2 <- plot.efficiency(distribution = distribution, method = method, 
                      size = c(100,250,500), plot = FALSE, etl=TRUE, 
                      replicates = 10000)

grid.arrange(p1, p2, nrow=2,ncol=1)

