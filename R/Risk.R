risk.calc <- function(beta, method = NA, distribution = NA, 
                      param = NA, data = NA, etl=FALSE) {
    
    if (is.na(method)) return(NA)
    else {
        if (method == "large-sample" && is.na(distribution)) return(NA)
        if (method == "large-sample" && distribution == "t" && is.na(param)) return(NA)
        if (method == "large-sample" && distribution == "skew-t" && is.na(param)) return(NA)
        if (method == "small-sample" && all(is.na(data))) return(NA)
    }
    
    norm.const <- 1; nu <- Inf; alpha <- 0
    ordinary <- modified <- c()
        
    moments <-         
        if(method == "large-sample") {            
            if (distribution == "gaussian") central.moments("gaussian") 
            else if (distribution == "t") {
                if ("nu" %in% names(param))  { 
                    nu <- param$nu 
                    norm.const <- sqrt((nu-2)/nu)
                    central.moments("t",list(nu=nu))                    
                }   else NA        
            } else if (distribution == "skew-t") {
                if (all(c("alpha","nu") %in% names(param)))  { 
                    alpha <-  param$alpha; nu <- param$nu 
                    norm.const <- 1
                    central.moments("skew-t",list(alpha=alpha, nu=nu))                    
                }   else NA        
            }
        }
    else if (method == "small-sample") 
        central.moments("small-sample", data=data)
    
    if (method == "small-sample" && distribution == "t") {
        mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
        fit <- fitdistr(data, mydt, start = list(df = 9), m = 0, s = 1, 
                        lower = 9, upper = 30)
        param$nu <- nu <- as.integer(fit$estimate["df"])                    
        norm.const <- sqrt((nu-2)/nu)
    }
    
    if (method == "small-sample" && distribution == "skew-t") {
        
        loglik.skewt <- function(x, param) {
            alpha <-  param[1]; nu <- param[2] 
            
            delta <- alpha/sqrt(1+alpha^2)
            b <- sqrt(nu/pi) * gamma(0.5*(nu - 1))/gamma(0.5*nu)
            
            term1 <- nu/(nu - 2); term2 <- b*delta 
            den <- sqrt(term1 - term2^2)
            xi <- -term2/den; omega <- 1/den
            
            -sum(log(dst(x, dp=c(xi=xi, omega=omega, alpha=alpha, nu=nu))))
        }
        
        start <- c(0, 9); lb <- c(-100, 9); ub <- c(100, 30)
        
        fit.skewt <- 
            optim(par = start, fn = loglik.skewt, x = data, method = "L-BFGS-B", 
                  lower = lb, upper = ub, control=list(maxit=1000))
        
        param$alpha <- alpha <- fit.skewt$par[1]
        param$nu    <- nu    <- fit.skewt$par[2]        
    }
    
    quant <- risk.quantile(beta, distribution, param)
    ordinary <- risk.ordinary(beta, method=method, distribution=distribution, 
                              param=list(alpha = alpha, nu=nu), quantile=quant, 
                              moments=moments, etl=etl, C=norm.const)
    
    quant <- risk.quantile(beta, "gaussian", param)
    modified <- risk.modified(beta,  method=method, distribution=distribution, 
                              quantile=quant, moments=moments, etl = etl)
    
    list(ordinary=ordinary, modified=modified, beta=beta, param=param, etl=etl)
}

risk.largesample.gaussian.efficiency <- function(etl=FALSE) {
    
    alphas  <- seq(0.01,0.05,by=0.04)
    main.grid <- c()
    
    for(alpha in alphas)
    {
        val <- risk.calc(alpha, method="large-sample", distribution="gaussian", etl=etl)
        bias.est <-  100*(val$ordinary$est - val$modified$est)/val$ordinary$est
        efficiency.se  <-  100*val$ordinary$se/val$modified$se
        grid <- matrix(c(alpha, bias.est, efficiency.se), nrow = 1)
        colnames(grid) <- c("sign", "est.bias", "se.efficiency")         
        main.grid <- rbind(main.grid,grid)
    }
    
    data.frame(main.grid)
}

risk.largesample.t.efficiency <- function(etl=FALSE) {
    
    alphas  <- seq(0.01,0.05,by=0.04)
    nus <- seq(from = 9, to = 30, by = 1)
    
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
                                   100*(val$ordinary$est - val$modified$est)/val$ordinary$est
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

risk.largesample.skewt.efficiency <- function(etl=FALSE) {
    
    betas  <- seq(0.01,0.05,by=0.04)
    nus <- seq(from = 9, to = 30, by = 1)
    
    main.grid <- c()
    st.alphas <- c(-20,10^-300,5,20)
    
    for(alpha in st.alphas) {
        for(beta in betas) {
            
            val <- sapply(nus,  
                          function(nu) {
                              param <- list(alpha = alpha, nu=nu)
                              val <- risk.calc(beta, method="large-sample", 
                                               distribution="skew-t", 
                                               param = param, etl = etl)
                              bias.est <-  
                                  100*(val$ordinary$est - val$modified$est)/val$ordinary$est
                              efficiency.se  <-  
                                  100*val$ordinary$se/val$modified$se
                              c(round(alpha), nu, beta, bias.est, efficiency.se,
                                val$ordinary$est,val$modified$est)
                          })
            
            grid <- t(val)        
            colnames(grid) <- c("slant", "df","sign", "est.bias", "se.efficiency", 
                                "VaR", "mVaR")         
            main.grid <- rbind(main.grid,grid)
        }
    }
    
    data.frame(main.grid)
}

risk.smallsample.gaussian.efficiency <- function(etl=FALSE, seed=1234, size) {
    
    set.seed(seed)
    norm.data <- sapply(size,rnorm)
    
    temp.func <- function(data) {
        
        alphas  <- seq(0.01,0.05,by=0.04)
        main.grid <- c()
        
        for(alpha in alphas)
        {
            val <- risk.calc(alpha, method="small-sample", 
                             distribution="gaussian", data=data, etl = etl)
            bias.est <-  100*(val$ordinary$est - val$modified$est)/val$ordinary$est
            efficiency.se  <-  100*val$ordinary$se/val$modified$se
            grid <- matrix(c(alpha, bias.est, efficiency.se), nrow = 1)
            colnames(grid) <- c("sign", "est.bias", "se.efficiency")         
            main.grid <- rbind(main.grid,grid)
        }
        
        main.grid
    }
    
    x <- lapply(norm.data, temp.func)
    main.grid <- lapply(1:length(size), function(i) cbind(size=size[i],x[[i]]))
    main.grid <- do.call(rbind,main.grid)    
    data.frame(main.grid)    
}

risk.smallsample.t.efficiency <- function(etl=FALSE, seed=99999, size) {

    temp.func <- function(size) {    
        
        alphas  <- seq(0.01,0.05,by=0.04)
        nus <- seq(from = 9, to = 30, by = 1)    
        main.grid <- c()
        
        for(alpha in alphas)
        {
            val <- sapply(nus,  
                          function(nu) {
                              set.seed(seed)    
                              data <- rt(size, df=nu)
                              val <- risk.calc(alpha, method="small-sample",
                                               distribution="t", 
                                               data=data, param=list(nu=nu),
                                               etl = etl)
                              bias.est <-  
                                  100*(val$ordinary$est - val$modified$est)/val$ordinary$est
                              efficiency.se  <-  
                                  100*val$ordinary$se/val$modified$se
                              c(size, nu, alpha, bias.est, efficiency.se)
                          })
            
            grid <- t(val)        
            colnames(grid) <- c("size", "df","sign", "est.bias", 
                                "se.efficiency")         
            main.grid <- rbind(main.grid,grid)
        }
        
        main.grid
    }
    
    main.grid <- lapply(size, temp.func)
    main.grid <- do.call(rbind,main.grid)    
    data.frame(main.grid)  
}

risk.smallsample.skewt.efficiency <- function(etl=FALSE, seed=99999, size) {
    
    temp.func <- function(size) {    
        betas  <- seq(0.01,0.05,by=0.04)
        nus <- seq(from = 9, to = 30, by = 1)    
        main.grid <- c()
        
        st.alphas <- c(-20,10^-300,5,20)
        
        for(alpha in st.alphas) {
            for(beta in betas)
            {
                val <- sapply(nus,  
                              function(nu) {
                                  set.seed(seed)    
                                                                                               
                                  delta <- alpha/sqrt(1+alpha^2)
                                  b <- sqrt(nu/pi) * 
                                      gamma(0.5*(nu - 1))/gamma(0.5*nu)
                                  
                                  term1 <- nu/(nu - 2); term2 <- b*delta 
                                  den <- sqrt(term1 - term2^2)
                                  xi <- -term2/den; omega <- 1/den
                                  
                                  data <- rst(size, xi = xi, omega = omega, 
                                              alpha = alpha, nu = nu)
                                  val <- risk.calc(beta, method="small-sample",
                                                   distribution="skew-t", 
                                                   data=data, 
                                                   param=list(alpha=alpha, nu=nu),
                                                   etl = etl)
                                  bias.est <-  
                                      100*(val$ordinary$est - val$modified$est)/val$ordinary$est
                                  efficiency.se  <-  
                                      100*val$ordinary$se/val$modified$se
                                  c(size, round(alpha), nu, beta, bias.est, 
                                    efficiency.se)
                              })
                
                grid <- t(val)        
                colnames(grid) <- c("size", "slant", "df","sign", "est.bias", 
                                    "se.efficiency")         
                main.grid <- rbind(main.grid,grid)
            }
        }
        
        main.grid
    }
    
    main.grid <- lapply(size, temp.func)
    main.grid <- do.call(rbind,main.grid)    
    data.frame(main.grid)  
}