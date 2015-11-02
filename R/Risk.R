risk.calc <- function(beta, method = NA, distribution = NA, 
                      param = NA, data = NA, etl=FALSE) {
    
    if (is.na(method)) return(NA)
    else {
        if (method == "large-sample" && is.na(distribution)) return(NA)
        if (method == "large-sample" && distribution == "t" && is.na(param)) return(NA)
        if (method == "small-sample" && all(is.na(data))) return(NA)
    }
    
    nu <- if ("nu" %in% names(param)) param$nu else Inf
    ordinary <- modified <- c()
        
    moments <-         
        if(method == "large-sample") {            
            if (distribution == "gaussian") central.moments("gaussian") 
            else if (distribution == "t") {
                if ("nu" %in% names(param))  { 
                    nu <- param$nu                     
                    central.moments("t",list(nu=nu))                    
                }   else NA        
            } 
        }
 
    if (method == "small-sample" && distribution == "gaussian") {
        moments <- c(mean(data), var(data))
        names(moments) <- c("mu", "mu2")
    } 
    
    if (method == "small-sample" && distribution == "t") {
        fit <- fitdistr(data, "t", df=nu)
        moments <- fit$estimate; moments[2] <- (nu/(nu-2))*(moments[2])^2
        names(moments) <- c("mu", "mu2")
    } 
        
    
    quant <- risk.quantile(beta, distribution, param)
    ordinary <- risk.ordinary(beta, method=method, distribution=distribution, 
                              param=list(nu=nu), quantile=quant, 
                              moments=moments, etl=etl)
    
    if (method == "small-sample") 
        moments <- central.moments("small-sample", data=data)
    
    quant <- risk.quantile(beta, "gaussian", param)
    modified <- risk.modified(beta,  method=method, distribution=distribution, 
                              quantile=quant, moments=moments, etl = etl)
    
    list(ordinary=ordinary, modified=modified, beta=beta, param=param, etl=etl)
}

risk.largesample.gaussian.efficiency <- function(etl=FALSE) {
    
    alphas  <- seq(0.01,0.05,by=0.001)
    main.grid <- c()
    
    for(alpha in alphas)
    {
        val <- risk.calc(alpha, method="large-sample", distribution="gaussian", etl=etl)
        bias.est <-  100*(val$modified$est - val$ordinary$est)/val$ordinary$est
        efficiency.se  <-  100*val$ordinary$se/val$modified$se
        grid <- matrix(c(alpha, bias.est, efficiency.se, val$ordinary$se,val$modified$se), nrow = 1)
        colnames(grid) <- c("sign", "est.bias", "se.efficiency", "ordinary", "modified")         
        main.grid <- rbind(main.grid,grid)
    }
    
    data.frame(round(main.grid,6))
}

risk.largesample.t.efficiency <- function(etl=FALSE) {
    
    alphas  <- c(0.01, 0.025, 0.05)
    nus <- seq(from = 9, to = 50, by = 1)  
    main.grid <- c()
    
    for(alpha in alphas)
    {
        val <- sapply(nus,  
                           function(nu) {
                               param <- list(nu=nu)
                               val <- risk.calc(alpha, method="large-sample", 
                                            distribution="t", param = param, 
                                            etl = etl)
                               bias.est <-  
                                   100*(val$modified$est - val$ordinary$est)/val$ordinary$est
                               efficiency.se  <-  
                                   100*val$ordinary$se/val$modified$se
                               c(nu, alpha, bias.est, efficiency.se)
                           })
        
        grid <- t(val)        
        colnames(grid) <- c("df","sign", "est.bias", "se.efficiency")         
        main.grid <- rbind(main.grid,grid)
    }
    
    data.frame(main.grid)
}



risk.smallsample.gaussian.efficiency <- function(etl=FALSE, seed=1234, replicates = 2500, size) {
    
    cl <- makeCluster(detectCores())
    clusterEvalQ(cl, library(foreach))
    clusterEvalQ(cl, library(moments))
    clusterSetRNGStream(cl, iseed = seed)
    
    registerDoSNOW(cl)
    
    alphas  <- seq(0.01,0.05,by=0.001)
    
    grid <- foreach(samplesize =  size, .combine = rbind ) %do% {
        
        M <- matrix(rnorm(samplesize*replicates), nrow = replicates)
        
        blocks <- foreach(alpha =  alphas, .combine = rbind ) %dopar% {
            
            #MLE estimator being unbiased we don't need small sample estimates
            val1 <- risk.calc(alpha, method="large-sample", 
                              distribution="gaussian", etl=etl)
            
            rows <- foreach(i =  1:replicates, .combine = rbind ) %do% {
                
                val2 <- risk.calc(alpha, method="small-sample", 
                                 distribution="gaussian", data=M[i,], etl = etl)
                
                c(val2$ordinary$est, val2$modified$est)
            }
            
            estOrdinary <- val1$ordinary$est; estModified <- rows[,2]
            bias <- mean(100*(estModified - estOrdinary)/estOrdinary)
            
            estOrdinary <- rows[,1]
            efficiency <- 100*sd(estOrdinary)/sd(estModified)
            c(bias, efficiency)
        }

        data.frame(size = rep(samplesize, length(alphas)), sign = alphas,
                   est.bias = blocks[,1], se.efficiency = blocks[,2])
    }
    
    stopCluster(cl)
    rownames(grid) <- NULL
    grid
}


risk.smallsample.t.efficiency <- function(etl=FALSE, seed=99999, size, 
                                          replicates = 2500) {
    
    set.seed(seed)
    
    alphas  <- c(0.01, 0.025, 0.05)
    
    nus <- seq(from = 9, to = 50, by = 1)  
    
    cl <- makeCluster(detectCores())
    clusterEvalQ(cl, library(foreach))
    clusterEvalQ(cl, library(moments))
    clusterEvalQ(cl, library(MASS))
    clusterSetRNGStream(cl, iseed = seed)
    
    registerDoSNOW(cl)
    
    temp.func <- function(nu) { 
        
        main.grid <- lapply(size, function(samplesize) {
            
            M <- matrix(rt(samplesize*replicates, df=nu), nrow = replicates)
            
            val <- foreach(alpha =  alphas, .combine = rbind ) %do% {
                
                val1 <- risk.calc(alpha, method="large-sample", 
                                  distribution="t", param = list(nu=nu), 
                                  etl = etl)
                rows <- 
                    foreach(i =  1:replicates, .combine = rbind) %do% {
                        
                        val2 <- risk.calc(alpha, method="small-sample",
                                         distribution="t", 
                                         data= M[i,], param=list(nu=nu),
                                         etl = etl)
                        
                        c(val2$ordinary$est, val2$modified$est)
                    }
                
                
                estOrdinary <- val1$ordinary$est; estModified <- rows[,2]
                bias <- mean(100*(estModified - estOrdinary)/estOrdinary)
                
                estOrdinary <- rows[,1]
                efficiency <- 100*sd(estOrdinary)/sd(estModified)
                c(samplesize, nu, alpha, bias, efficiency)
            }
            
            rownames(val) <- NULL
            colnames(val) <- c("size", "df","sign", "est.bias", 
                               "se.efficiency")         
            val
        })
        
        main.grid <- do.call(rbind, main.grid)
        rownames(main.grid) <- NULL
        main.grid
    }
    
    main.grid <- parLapply(cl, nus, temp.func)
    stopCluster(cl)
    
    main.grid <- do.call(rbind,main.grid)    
    data.frame(main.grid)  
}