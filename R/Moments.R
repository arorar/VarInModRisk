# assume mean=0; variance=1
gaussian.central.moment <- function(order) {
    
    if(order %in% c(1,3,5,7)) 0
    else if (order == 2) 1
    else if (order == 4) 3
    else if (order == 6) 15
    else if (order == 8) 105
}

# assume mean=0; variance=1
student.t.central.moment <- function(order, param) {
    
    if ("nu" %in% names(param))  { 
        nu <- param$nu 
        
        if(order %in% c(1,3,5,7)) 0
        else if (order == 2) nu/(nu - 2)
        else if (order == 4) (3*nu^2)/((nu - 2)*(nu - 4))
        else if (order == 6) (15*nu^3)/((nu - 2)*(nu - 4)*(nu - 6))
        else if (order == 8) (105*nu^4)/((nu - 2)*(nu - 4)*(nu - 6)*(nu - 8))
    } else NA
}

empirical.central.moment <- function(R, order) {
    
    moment(R, order=order, central=TRUE)    
}

central.moments <- function(type=c("gaussian","t","small-sample"), 
                            param = NA, data = NA) {
    
    if (type == "t" && is.na(param)) return(NA)
    else if (type == "small-sample" && is.na(data)) return(NA)
    
    temp.moments <- function(order) {
        
        if (!type %in% c("gaussian","t","small-sample")) NA
        
        if (type == "gaussian") gaussian.central.moment(order)
        else if (type == "t")   student.t.central.moment(order, param) 
        else if (type == "small-sample") empirical.central.moment(data, order)
        else NA
    }
    
    moments <- sapply(1:8, temp.moments)
    names(moments) <- c("mu", paste("mu", 2:8, sep=""))
    moments    
}