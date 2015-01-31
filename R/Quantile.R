risk.quantile <- function(alpha, distribution = NA, param = NA) {
    
    
    
    if (is.na(distribution)) NA
    if (distribution == "t" && is.na(param)) NA        
    
    if (distribution == "gaussian") qnorm(alpha)
    else if (distribution == "t") {            
        if ("nu" %in% names(param))  { 
            nu <- param$nu 
            qt(p = alpha, df = nu)
        }   else NA
    }
}