risk.quantile <- function(beta, distribution = NA, param = NA) {    
    
    if (is.na(distribution)) NA
    if (distribution == "t" && is.na(param)) NA        
    
    if (distribution == "gaussian") qnorm(beta)
    else if (distribution == "t") {            
        if ("nu" %in% names(param))  { 
            nu <- param$nu             
            qt(p = beta, df = nu)
        }   else NA
    } 
}