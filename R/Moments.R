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

skew.t.central.moment <- function(order, param) {
        
    if (!all(c("alpha","nu") %in% names(param)))  return(NA)
    
    alpha <-  param$alpha;  nu <- param$nu 
    
    delta <- alpha/sqrt(1+alpha^2)
    b <- sqrt(nu/pi) * gamma(0.5*(nu - 1))/gamma(0.5*nu)
    
    term1 <- nu/(nu - 2); term2 <- b*delta 
    den <- sqrt(term1 - term2^2)
    xi <- -term2/den; omega <- 1/den
        
    if(order == 1) {
        xi + b*delta*omega
    }
    else if(order == 2) {
        (omega^2*(nu + 2*b^2*delta^2 - b^2*delta^2*nu))/(nu - 2)
    }
    else if(order == 3) {
        (b*(2*b^2*nu^2 - 10*b^2*nu + 12*b^2 - nu^2 + 2*nu)*delta^3*omega^3 + 
             3*b*nu*delta*omega^3)/((nu - 2)*(nu - 3))
    } 
    else if(order == 4) {
        (3*nu^2*omega^4*(nu - 3) - 6*b^2*delta^2*nu*omega^4*(nu^2 - 
          5*nu + 4) + b^2*delta^4*omega^4*(4*nu - 3*b^2*nu + 9*b^2)*
             (nu^2 - 6*nu + 8))/((nu - 2)*(nu - 3)*(nu - 4))
    } 
    else if(order == 5) {
        (15*b*delta*nu^2*omega^5*(2*nu - 7) + b*delta^5*omega^5*
             (nu^2 - 6*nu + 8)*(4*b^4*nu^2 - 32*b^4*nu + 60*b^4 - 10*b^2*nu^2 + 
            50*b^2*nu + 3*nu^2) + 10*b*delta^3*nu*omega^5*(nu - 4)*(2*b^2*nu^2 - 
            13*b^2*nu + 15*b^2 - nu^2 + 2*nu))/((nu - 2)*(nu - 3)*(nu - 4)*(nu - 5))
    }
    else if(order == 6) {
        -(b^2*delta^6*omega^6*(nu^3 - 12*nu^2 + 44*nu - 48)*(5*b^4*nu^2 - 
        40*b^4*nu + 75*b^4 - 20*b^2*nu^2 + 100*b^2*nu + 18*nu^2) - 
        15*nu^3*omega^6*(nu^2 - 8*nu + 15) + 
            45*b^2*delta^2*nu^2*omega^6*(nu^3 - 10*nu^2 + 25*nu - 6) + 
            15*b^2*delta^4*nu*omega^6*(nu^2 - 10*nu + 24)*
            (3*b^2*nu^2 - 20*b^2*nu + 25*b^2 - 4*nu^2 + 8*nu))/
            ((nu - 2)*(nu - 3)*(nu - 4)*(nu - 5)*(nu - 6))
    }
    else if(order == 7) {
        -(b*delta^7*omega^7*(nu^3 - 12*nu^2 + 44*nu - 48)*
              (- 6*b^6*nu^3 + 90*b^6*nu^2 - 426*b^6*nu + 630*b^6 + 35*b^4*nu^3 - 
        420*b^4*nu^2 + 1225*b^4*nu - 63*b^2*nu^3 + 441*b^2*nu^2 + 15*nu^3) - 
            315*b*delta*nu^3*omega^7*(nu^2 - 9*nu + 19) + 
            21*b*delta^5*nu*omega^7*(nu^2 - 10*nu + 24)*
            (- 4*b^4*nu^3 + 55*b^4*nu^2 - 224*b^4*nu + 245*b^4 + 10*b^2*nu^3 - 
                 90*b^2*nu^2 + 140*b^2*nu - 3*nu^3 + 6*nu^2) + 
            105*b*delta^3*nu^2*omega^7*(nu - 6)*(- 2*b^2*nu^3 + 24*b^2*nu^2 - 
                79*b^2*nu + 63*b^2 + nu^3 - 6*nu^2 + 8*nu))/
            ((nu - 2)*(nu - 3)*(nu - 4)*(nu - 5)*(nu - 6)*(nu - 7))
    }
    else if(order == 8) {
        (105*nu^4*omega^8*(nu^3 - 15*nu^2 + 71*nu - 105) + 
             420*b^2*delta^2*nu^3*omega^8*(- nu^4 + 17*nu^3 - 89*nu^2 + 127*nu + 72) +
             b^2*delta^8*omega^8*(nu^4 - 20*nu^3 + 140*nu^2 - 400*nu + 384)*
             (- 7*b^6*nu^3 + 105*b^6*nu^2 - 497*b^6*nu + 735*b^6 + 56*b^4*nu^3 - 
                  672*b^4*nu^2 + 1960*b^4*nu - 168*b^2*nu^3 + 1176*b^2*nu^2 + 
            120*nu^3) + 210*b^2*delta^4*nu^2*omega^8*(nu^2 - 14*nu + 48)*
             (- 3*b^2*nu^3 + 37*b^2*nu^2 - 129*b^2*nu + 119*b^2 + 4*nu^3 - 
                  24*nu^2 + 32*nu) + 28*b^2*delta^6*nu*omega^8*(nu^3 - 18*nu^2 + 
        104*nu - 192)*(- 5*b^4*nu^3 + 69*b^4*nu^2 - 283*b^4*nu + 315*b^4 + 
            20*b^2*nu^3 - 180*b^2*nu^2 + 280*b^2*nu - 18*nu^3 + 36*nu^2))/
            ((nu - 2)*(nu - 3)*(nu - 4)*(nu - 5)*(nu - 6)*(nu - 7)*(nu - 8))
    }         
}

empirical.central.moment <- function(R, order) {
    
    moment(R, order=order, central=TRUE)    
}

central.moments <- function(type=c("gaussian","t","skew-t","small-sample"), 
                            param = NA, data = NA) {
    
    if (type == "t" && is.na(param)) return(NA)
    else if (type == "small-sample" && is.na(data)) return(NA)
    
    temp.moments <- function(order) {
        
        if (!type %in% c("gaussian","t","skew-t","small-sample")) NA
        
        if (type == "gaussian") gaussian.central.moment(order)
        else if (type == "t")   student.t.central.moment(order, param)
        else if (type == "skew-t")   skew.t.central.moment(order, param)
        else if (type == "small-sample") empirical.central.moment(data, order)
        else NA
    }
    
    moments <- sapply(1:8, temp.moments)
    names(moments) <- c("mu", paste("mu", 2:8, sep=""))
    moments    
}