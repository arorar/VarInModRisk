risk.factor <- function(z, distribution, moments, param = NA) {
    if(distribution == "gaussian") dnorm(z)
    else if(distribution == "t") {
        if ("nu" %in% names(param))  { 
            nu <- param$nu            
            dt(z, df = nu - 2)
        } else NA
    } else if(distribution == "skew-t") {
        if (all(c("alpha","nu") %in% names(param)))  { 
            alpha <-  param$alpha; nu <- param$nu
            
            delta <- alpha/sqrt(1+alpha^2)
            b <- sqrt(nu/pi) * gamma(0.5*(nu - 1))/gamma(0.5*nu)
            
            term1 <- nu/(nu - 2); term2 <- b*delta 
            den <- sqrt(term1 - term2^2)
            xi <- -term2/den; omega <- 1/den
            
            mu  <- moments["mu"];  sigma <- sqrt(moments["mu2"])            
            var <- mu + sigma*z
            
            integrand <- function(x, xi, omega, alpha, nu) {
                x*dst(x,dp = c(xi = xi, omega = omega, alpha = alpha, nu = nu))
            }
            
            -integrate(integrand, lower = -Inf, upper = var, xi = xi, omega = omega, 
                      alpha = alpha, nu =nu)$value
        } else NA
    }    
}

risk.ordinary <- function(beta, method, distribution, param, quantile, moments, 
                          etl) {    

    z   <- quantile; dz <- risk.factor(z, distribution, moments, param)
    mu  <- moments["mu"];  mu2 <- moments["mu2"]
    mu3 <- moments["mu3"]; mu4 <- moments["mu4"]
    
    sigma <- sqrt(mu2); skew <- mu3/mu2^(3/2); exkurt <- mu4/mu2^2 - 3
    
    est <- if(!etl) mu + sigma*z else mu - sigma/beta*dz        
    
    variance <- 
        if(!etl) 
            (z^2/2*(1+exkurt/2) + skew*z + 1 )*sigma^2 
        else
            ((dz^2*(exkurt + 2))/(4*beta^2) - (dz*skew)/beta + 1)*sigma^2 
    
    names(est) <- names(variance) <- NULL
    
    list(est = est, se = sqrt(variance))
}

risk.modified <- function(beta, method, distribution, quantile, moments, etl) {

    z   <-  quantile; dz <- dnorm(z)
    mu  <-  moments["mu"];  mu2 <- moments["mu2"]
    mu3 <-  moments["mu3"]; mu4 <- moments["mu4"]
    mu5 <-  moments["mu5"]; mu6 <- moments["mu6"]
    mu7 <-  moments["mu7"]; mu8 <- moments["mu8"]
        
    sigma <- sqrt(mu2) 
    skew <- mu3/mu2^(3/2) 
    exkurt <- mu4/mu2^2 - 3
    g <- z + (z^2-1)*skew/6 + (z^3-3*z)*exkurt/24 - (2*z^3-5*z)*skew^2/36
    dg <- dnorm(g)
    
    est <- if (!etl) {
        mvar <- mu + sigma*g
        mvar
    } else {
        term1 <-  (g^3)
        term2 <-  (g^4 - 2*g^2 - 1)
        term3 <-  (g^6 - 9*g^4 + 9*g^2 + 3)        
        
        factor <- 1 + (term1)*skew/6 + (term2)*exkurt/24 + (term3)*skew^2/72
        metl <- mu - dg*sigma/beta*factor
        metl
    }
        
    variance <-  if(!etl) 
        -(sigma^8*((- 81*exkurt^3 - 882*exkurt^2 - 3132*exkurt - 3564)*z^6 + 
        (486*exkurt^3 + 6588*exkurt^2 + 30888*exkurt + 39528)*z^4 + 
        (- 729*exkurt^3 - 11826*exkurt^2 - 73116*exkurt - 112428)*z^2 +
        10368*exkurt + 5184) - 
        skew*(144*mu5*sigma^3*(15*exkurt*z^2 - 11*exkurt*z^4 + 
        2*exkurt*z^6 + 128*z^2 - 78*z^4 + 10*z^6 - 8) - 
        96*sigma*z*(15*mu7*z - 11*mu7*z^3 + 2*mu7*z^5 + 29*mu6*sigma - 
        40*mu6*z^2*sigma + 11*mu6*z^4*sigma) + 
        144*sigma^7*z*sigma*(3*exkurt^2*z^4 - 12*exkurt^2*z^2 + 
        9*exkurt^2 + 54*exkurt*z^4 - 274*exkurt*z^2 + 304*exkurt + 
        120*z^4 - 768*z^2 + 1080)) - 324*mu8*z^2 + 216*mu8*z^4 - 
        36*mu8*z^6 + skew^3*(160*mu5*sigma^3*z^2*(2*z^2 - 5)^2 + 
        96*sigma^7*z*sigma*(25*exkurt - 35*exkurt*z^2 + 10*exkurt*z^4 - 
        306*z^2 + 72*z^4 + 324)) - skew^2*(sigma^8*((- 360*exkurt^2 - 
        3696*exkurt - 7632)*z^6 + (1980*exkurt^2 + 22800*exkurt + 
        56376)*z^4 + (- 2700*exkurt^2 - 33504*exkurt - 95256)*z^2 + 
        576*exkurt + 12096) + 4320*mu5*z*sigma^3 - 6048*mu5*z^3*sigma^3 + 
        1728*mu5*z^5*sigma^3 + 8*mu6*sigma^2*z^2*(62*z^4 - 325*z^2 + 425)) + 
        144*z^5*sigma*(18*mu5*sigma^2 - 2*mu7 + 3*exkurt*mu5*sigma^2) - 
        576*z^3*sigma*(27*mu5*sigma^2 - 2*mu7 + 3*exkurt*mu5*sigma^2) + 
        36*mu6*sigma^2*(27*exkurt*z^2 - 18*exkurt*z^4 + 3*exkurt*z^6 + 212*z^2 -
        112*z^4 + 12*z^6 - 16) + 432*z*sigma*(38*mu5*sigma^2 - 2*mu7 + 
        3*exkurt*mu5*sigma^2) + 4*skew^4*sigma^8*z^2*(2*z^2 - 5)*(125*exkurt -
        50*exkurt*z^2 - 268*z^2 + 610))/(20736*sigma^6)
    else
        (36*dg^2*mu8 + 20736*beta^2*sigma^8 + 26892*dg^2*sigma^8 + 
             16956*dg^2*exkurt*sigma^8 - 1296*dg^2*mu6*sigma^2 + 
             144*dg^2*mu8*g^2 + 72*dg^2*mu8*g^4 - 144*dg^2*mu8*g^6 + 
             36*dg^2*mu8*g^8 - 12492*dg^2*skew^2*sigma^8 + 
             1026*dg^2*skew^4*sigma^8 + 2178*dg^2*exkurt^2*sigma^8 + 
             81*dg^2*exkurt^3*sigma^8 + 40176*dg^2*sigma^8*g^2 - 
             5832*dg^2*sigma^8*g^4 - 19440*dg^2*sigma^8*g^6 + 
             3564*dg^2*sigma^8*g^8 - 4536*dg^2*skew^2*exkurt*sigma^8 + 
             225*dg^2*skew^4*exkurt*sigma^8 + 324*dg^2*skew^2*mu6*sigma^2 - 
             360*dg^2*skew^3*mu5*sigma^3 + 29808*dg^2*exkurt*sigma^8*g^2 - 
             2376*dg^2*exkurt*sigma^8*g^4 - 15984*dg^2*exkurt*sigma^8*g^6 + 
             3132*dg^2*exkurt*sigma^8*g^8 - 3456*dg^2*mu6*sigma^2*g^2 + 
             2304*dg^2*mu6*sigma^2*g^6 - 432*dg^2*mu6*sigma^2*g^8 - 
             144*dg^2*skew*mu7*sigma - 270*dg^2*skew^2*exkurt^2*sigma^8 - 
             53964*dg^2*skew^2*sigma^8*g^2 + 19512*dg^2*skew^2*sigma^8*g^4 + 
             6876*dg^2*skew^4*sigma^8*g^2 + 64296*dg^2*skew^2*sigma^8*g^6 + 
             3798*dg^2*skew^4*sigma^8*g^4 - 32940*dg^2*skew^2*sigma^8*g^8 - 
             24024*dg^2*skew^4*sigma^8*g^6 + 4500*dg^2*skew^2*sigma^8*g^10 + 
             16566*dg^2*skew^4*sigma^8*g^8 - 144*dg^2*skew^2*sigma^8*g^12 - 
             3252*dg^2*skew^4*sigma^8*g^10 + 194*dg^2*skew^4*sigma^8*g^12 + 
             6120*dg^2*exkurt^2*sigma^8*g^2 + 324*dg^2*exkurt^3*sigma^8*g^2 + 
             468*dg^2*exkurt^2*sigma^8*g^4 + 162*dg^2*exkurt^3*sigma^8*g^4 - 
             3528*dg^2*exkurt^2*sigma^8*g^6 - 324*dg^2*exkurt^3*sigma^8*g^6 + 
             882*dg^2*exkurt^2*sigma^8*g^8 + 81*dg^2*exkurt^3*sigma^8*g^8 + 
             6048*dg^2*mu5*g^3*sigma^3 - 288*dg^2*mu7*g^3*sigma + 
             5184*dg^2*mu5*g^5*sigma^3 - 576*dg^2*mu7*g^5*sigma - 
             2592*dg^2*mu5*g^7*sigma^3 + 288*dg^2*mu7*g^7*sigma + 
             1728*beta*dg*mu5*sigma^3 + 2736*dg^2*skew*mu5*sigma^3 - 
             108*dg^2*exkurt*mu6*sigma^2 - 20736*beta*dg*sigma^8*g^3 - 
             19656*dg^2*skew^2*exkurt*sigma^8*g^2 + 864*dg^2*skew^2*exkurt*sigma^8*g^4 + 
             1350*dg^2*skew^4*exkurt*sigma^8*g^2 + 29952*dg^2*skew^2*exkurt*sigma^8*g^6 + 
             675*dg^2*skew^4*exkurt*sigma^8*g^4 - 16200*dg^2*skew^2*exkurt*sigma^8*g^8 - 
             3900*dg^2*skew^4*exkurt*sigma^8*g^6 + 2376*dg^2*skew^2*exkurt*sigma^8*g^10 + 
             2475*dg^2*skew^4*exkurt*sigma^8*g^8 - 96*dg^2*skew^2*exkurt*sigma^8*g^12 - 
             450*dg^2*skew^4*exkurt*sigma^8*g^10 + 25*dg^2*skew^4*exkurt*sigma^8*g^12 + 
             1764*dg^2*skew^2*mu6*sigma^2*g^2 - 2160*dg^2*skew^3*mu5*sigma^3*g^2 + 
             792*dg^2*skew^2*mu6*sigma^2*g^4 - 1080*dg^2*skew^3*mu5*sigma^3*g^4 - 
             4056*dg^2*skew^2*mu6*sigma^2*g^6 + 6240*dg^2*skew^3*mu5*sigma^3*g^6 + 
             2244*dg^2*skew^2*mu6*sigma^2*g^8 - 3960*dg^2*skew^3*mu5*sigma^3*g^8 - 
             348*dg^2*skew^2*mu6*sigma^2*g^10 + 720*dg^2*skew^3*mu5*sigma^3*g^10 + 
             16*dg^2*skew^2*mu6*sigma^2*g^12 - 40*dg^2*skew^3*mu5*sigma^3*g^12 + 
             432*dg^2*exkurt*mu5*g^3*sigma^3 + 864*dg^2*exkurt*mu5*g^5*sigma^3 - 
             432*dg^2*exkurt*mu5*g^7*sigma^3 + 216*dg^2*skew*exkurt*mu5*sigma^3 -
             6912*beta*dg*exkurt*sigma^8*g^3 - 720*dg^2*skew*mu7*sigma*g^2 - 
             288*dg^2*skew*mu7*sigma*g^4 + 1248*dg^2*skew*mu7*sigma*g^6 - 
             528*dg^2*skew*mu7*sigma*g^8 + 48*dg^2*skew*mu7*sigma*g^10 - 
             1350*dg^2*skew^2*exkurt^2*sigma^8*g^2 - 540*dg^2*skew^2*exkurt^2*sigma^8*g^4 + 
             2340*dg^2*skew^2*exkurt^2*sigma^8*g^6 - 990*dg^2*skew^2*exkurt^2*sigma^8*g^8 + 
             90*dg^2*skew^2*exkurt^2*sigma^8*g^10 - 
             1296*dg^2*skew^2*mu5*g^3*sigma^3 - 3888*dg^2*skew^2*mu5*g^5*sigma^3 + 
             3888*dg^2*skew^2*mu5*g^7*sigma^3 - 432*dg^2*skew^2*mu5*g^9*sigma^3 - 
             36288*dg^2*skew*sigma^7*g^3*sigma - 36288*dg^2*skew*sigma^7*g^5*sigma +
             25920*dg^2*skew*sigma^7*g^7*sigma - 1728*dg^2*skew*sigma^7*g^9*sigma - 
             41472*beta*dg*skew*sigma^7*sigma + 6912*beta*dg*skew^2*sigma^8*g^3 + 
             3456*beta*dg*mu5*g^2*sigma^3 - 1728*beta*dg*mu5*g^4*sigma^3 + 
             10512*dg^2*skew*mu5*sigma^3*g^2 - 3168*dg^2*skew*mu5*sigma^3*g^4 - 
             10656*dg^2*skew*mu5*sigma^3*g^6 + 4464*dg^2*skew*mu5*sigma^3*g^8 - 
             432*dg^2*skew*mu5*sigma^3*g^10 - 432*dg^2*exkurt*mu6*sigma^2*g^2 - 
             216*dg^2*exkurt*mu6*sigma^2*g^4 + 432*dg^2*exkurt*mu6*sigma^2*g^6 - 
             108*dg^2*exkurt*mu6*sigma^2*g^8 + 4896*dg^2*skew^3*sigma^7*g^3*sigma + 
             15840*dg^2*skew^3*sigma^7*g^5*sigma - 16992*dg^2*skew^3*sigma^7*g^7*sigma + 
             2016*dg^2*skew^3*sigma^7*g^9*sigma + 4320*beta*dg*skew^3*sigma^7*sigma + 
             1080*dg^2*skew*exkurt*mu5*sigma^3*g^2 + 432*dg^2*skew*exkurt*mu5*sigma^3*g^4 - 
             1872*dg^2*skew*exkurt*mu5*sigma^3*g^6 + 792*dg^2*skew*exkurt*mu5*sigma^3*g^8 -
             72*dg^2*skew*exkurt*mu5*sigma^3*g^10 - 432*dg^2*skew*exkurt^2*sigma^7*g^3*sigma + 
             720*dg^2*skew^3*exkurt*sigma^7*g^3*sigma - 864*dg^2*skew*exkurt^2*sigma^7*g^5*sigma + 
             2160*dg^2*skew^3*exkurt*sigma^7*g^5*sigma + 432*dg^2*skew*exkurt^2*sigma^7*g^7*sigma - 
             2160*dg^2*skew^3*exkurt*sigma^7*g^7*sigma + 240*dg^2*skew^3*exkurt*sigma^7*g^9*sigma -
             51840*beta*dg*skew*sigma^7*g^2*sigma + 41472*beta*dg*skew*sigma^7*g^4*sigma - 
             3456*beta*dg*skew*sigma^7*g^6*sigma + 864*dg^2*skew*mu6*sigma*g^3*sigma + 
             2304*dg^2*skew*mu6*sigma*g^5*sigma - 2016*dg^2*skew*mu6*sigma*g^7*sigma + 
             192*dg^2*skew*mu6*sigma*g^9*sigma + 12960*beta*dg*skew^3*sigma^7*g^2*sigma - 
             12960*beta*dg*skew^3*sigma^7*g^4*sigma + 1440*beta*dg*skew^3*sigma^7*g^6*sigma - 
             10080*dg^2*skew*exkurt*sigma^7*g^3*sigma - 16704*dg^2*skew*exkurt*sigma^7*g^5*sigma +
             13536*dg^2*skew*exkurt*sigma^7*g^7*sigma - 1152*dg^2*skew*exkurt*sigma^7*g^9*sigma - 
             6048*beta*dg*skew*exkurt*sigma^7*sigma - 15552*beta*dg*skew*exkurt*sigma^7*g^2*sigma + 
             12960*beta*dg*skew*exkurt*sigma^7*g^4*sigma - 
             1152*beta*dg*skew*exkurt*sigma^7*g^6*sigma)/(20736*beta^2*sigma^6)
    
        names(est) <- names(variance) <- NULL
        list(est = est, se = sqrt(variance))
}