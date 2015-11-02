risk.factor <- function(z, distribution, moments, param = NA) {
    if(distribution == "gaussian") dnorm(z)
    else if(distribution == "t") {
        if ("nu" %in% names(param))  { 
            nu <- param$nu 
            lambda <- sqrt(nu/(nu-2))
            lambda*dt(z/lambda, df = nu - 2)
        } else NA
    } 
}

risk.ordinary <- function(beta, method, distribution, param, quantile, moments, 
                          etl) {    

    z   <- quantile; dz <- risk.factor(z, distribution, moments, param)
    mu  <- moments["mu"];  mu2 <- moments["mu2"]; sigma <- sqrt(mu2)
    
    if(distribution == "t") {
        nu <- param$nu
        sigma <- sqrt((nu-2)/nu)*sigma
    }
    
    est <- if(!etl) mu + sigma*z else mu - sigma/beta*dz        
    
    variance <-
        if(distribution == "gaussian") {
            if(!etl) sigma^2*(1+z^2/2) else sigma^2*(1+dz^2/(2*beta^2)) 
        }
        else if(distribution == "t") {
            nu <- param$nu
            if(!etl) sigma^2*((nu+3)/(nu+1) + (nu+3)/nu*z^2/2)
            else sigma^2*((nu+3)/(nu+1) + (nu+3)/nu*dz^2/(2*beta^2))
        }

    names(est) <- names(variance) <- NULL
    
    list(est = est, se = sqrt(variance))
}

risk.modified <- function(beta, method, distribution, quantile, moments, etl) {

    z   <-  quantile; dz <- dnorm(z)
    mu  <-  moments["mu"];  mu2 <- moments["mu2"]; sigma <- sqrt(mu2)
    mu3 <-  moments["mu3"]; mu4 <- moments["mu4"]; skew <- mu3/sigma^3; exkurt <- mu4/sigma^4 - 3
    mu5 <-  moments["mu5"]; mu6 <- moments["mu6"]
    mu7 <-  moments["mu7"]; mu8 <- moments["mu8"]
   
    V <- matrix(nrow = 4, ncol = 4)     
     
    V[1,] <- c(sigma^2,            mu3,                         mu4- 3*sigma^4,                              mu5 -4*mu3*sigma^2                                                 )
    V[2,] <- c(mu3,                mu4-sigma^4,                 mu5-4*sigma^2*mu3,                           mu6-sigma^2*mu4-4*mu3^2                              )
    V[3,] <- c(mu4- 3*sigma^4,     mu5-4*sigma^2*mu3,           mu6-mu3^2-6*sigma^2*mu4+9*sigma^6,           mu7-mu3*mu4-3*sigma^2*mu5-4*mu3*mu4+12*sigma^4*mu3   )
    V[4,] <- c(mu5 -4*mu3*sigma^2, mu6-sigma^2*mu4-4*mu3^2,     mu7-3*sigma^2*mu5-5*mu3*mu4+12*sigma^4*mu3,  mu8-mu4^2-8*mu3*mu5+16*sigma^2*mu3^2                 )
    
    C0 <- z; C1 <- (z^2 - 1)/6 ; C2 <- (z^3 - 3*z)/24; C3 <- -(2*z^3 - 5*z)/36
    g <- C0 + C1*skew + C2*exkurt + C3*skew^2
    dg <- dnorm(g)
    
    if (etl) {
        C0 <- -dg/beta; C1 <- -dg*g^3/(6*beta); C2 <- -dg*(g^4 - 2*g^2 - 1)/(24*beta)
        C3 <- -dg*(g^6 - 9*g^4 + 9*g^2 + 3)/(72*beta)
    }

    est <- mu + sigma*(C0 + C1*skew + C2*exkurt + C3*skew^2)

    D <- matrix(c( 1,0,0,0,0, 
                   0,1/(2*sigma), -mu3/sigma^4, -(3*sigma^4 + 3*mu4)/(2*sigma^5), -(5*mu3^2)/(2*sigma^7), 
                   0,0, 1/sigma^2,0,(2*mu3)/sigma^5, 
                   0,0, 0, 1/sigma^3,0), nrow=4, byrow = TRUE)
    
    x <- as.matrix(c(1, C0, C1, C2, C3))
    grad <- D %*% x 
    
    if (etl) {
        
        Dg <- as.matrix(c(0, 
                mu3^2*z*(2*z^2-5)/(12*sigma^8) - mu3*(z^2-1)/(4*sigma^5) - mu4*z*(z^2-3)/(12*sigma^6),
                (z^2-1)/(6*sigma^3) - mu3*z*(2*z^2-5)/(18*sigma^6),
                z*(z^2-3)/(24*sigma^4)))
        
        Dx <- as.matrix(c(0,
                g,
                1/6*(g^4 - 3*g^2),
                1/72*(g^7 - 15*g^5 + 45*g^3 - 15*g),
                1/24*(g^5 - 6*g^3 + 3*g)))
        
        y <- as.matrix(c(mu, sigma, sigma*skew, sigma*skew^2, sigma*exkurt))
        grad <- grad + dg*sigma/beta*Dg %*% t(Dx) %*% y
    }
    
    variance <-  t(grad) %*% V %*% grad

    names(est) <- names(variance) <- NULL
    list(est = est, se = sqrt(variance))
}