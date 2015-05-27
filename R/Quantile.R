risk.quantile <- function(beta, distribution = NA, param = NA) {    
    
    if (is.na(distribution)) NA
    if (distribution == "t" && is.na(param)) NA        
    if (distribution == "skew-t" && is.na(param)) NA        
    
    if (distribution == "gaussian") qnorm(beta)
    else if (distribution == "t") {            
        if ("nu" %in% names(param))  { 
            nu <- param$nu             
            qstd(p = beta, nu = nu)
        }   else NA
    } else if (distribution == "skew-t") {            
        if (all(c("alpha","nu") %in% names(param)))  { 
            alpha <-  param$alpha
            nu <- param$nu 
            
            delta <- alpha/sqrt(1+alpha^2)
            b <- sqrt(nu/pi) * gamma(0.5*(nu - 1))/gamma(0.5*nu)
            
            term1 <- nu/(nu - 2); term2 <- b*delta 
            den <- sqrt(term1 - term2^2)
            xi <- -term2/den; omega <- 1/den
            
            qst(p = beta, xi = xi, omega = omega, alpha = alpha, nu = nu)            
        }   else NA
    }
}