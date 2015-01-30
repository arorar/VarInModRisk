risk.quantile <- function(alpha, method = NA, distribution = NA, 
                          param = NA, data = NA) {
    
    if (is.na(method)) return(NA)
    else {
        if (method == "large-sample" && is.na(distribution)) NA
        if (method == "large-sample" && distribution == "t" && is.na(param)) NA
        if (method == "small-sample" && all(is.na(data))) NA
    }
    
    if (method == "large-sample" && distribution == "gaussian") qnorm(alpha)
    else if (method == "large-sample" && distribution == "t") {            
        if ("nu" %in% names(param))  { 
            nu <- param$nu 
            qt(p = alpha, df = nu)
        }   else NA
    } else if (method == "small-sample") { 
        x <- quantile(data, alpha)
        names(x) <- NULL
        x
    }
}