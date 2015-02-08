risk.factor <- function(z, distribution, moments, param = NA, C = NA) {
    if(distribution == "gaussian") dnorm(z)
    else if(distribution == "t") {
        if ("nu" %in% names(param))  { 
            nu <- param$nu
            dt(sqrt(C)*z, df = nu - 2)/sqrt(C)
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
            var <- mu + C*sigma*z
            
            integrand <- function(x, xi, omega, alpha, nu) {
                x*dst(x,dp = c(xi = xi, omega = omega, alpha = alpha, nu = nu))
            }
            
            -integrate(integrand, lower = -Inf, upper = var, xi = xi, omega = omega, 
                      alpha = alpha, nu =nu)$value
        } else NA
    }    
}

risk.ordinary <- function(beta, method, distribution, param, quantile, moments, 
                          etl, C) {    

    z   <- quantile; dz <- risk.factor(z, distribution, moments, param, C = C)
    mu  <- moments["mu"];  mu2 <- moments["mu2"]
    mu3 <- moments["mu3"]; mu4 <- moments["mu4"]
    
    sigma <- sqrt(mu2); skew <- mu3/mu2^(3/2); exkurt <- mu4/mu2^2 - 3
    
    est <- if(!etl) mu + C*sigma*z else mu - C*sigma/beta*dz        
    
    variance <- 
        if(!etl) 
            ((C*z)^2/2*(1+exkurt/2) + skew*C*z + 1 )*sigma^2 
        else
            (((C*dz)^2*(exkurt + 2))/(4*beta^2) - (C*dz*skew)/beta + 1)*sigma^2 
    
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
    
    est <- if (!etl) {
        factor <- z + (z^2-1)*skew/6 + (z^3-3*z)*exkurt/24 - (2*z^3-5*z)*skew^2/36
        mvar <- mu + sigma*factor
        mvar
    } else {
        term1 <-  (dz*z^3)
        term2 <-  (z^4 - 2*z^2 - 1)*dz
        term3 <-  (z^6 - 9*z^4 + 9*z^2 + 3)*dz        
        
        factor <- dz + (term1)*skew/6 + (term2)*exkurt/24 + (term3)*skew^2/72
        metl <- mu - sigma/beta*factor
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
        (20736*beta^2*sigma^8 + 1440*beta*dz*skew^3*sigma^8*z^6 - 
             12960*beta*dz*skew^3*sigma^8*z^4 + 12960*beta*dz*skew^3*sigma^8*z^2 + 
             4320*beta*dz*skew^3*sigma^8 + 6912*beta*dz*skew^2*sigma^8*z^3 - 
             1152*beta*dz*skew*exkurt*sigma^8*z^6 + 12960*beta*dz*skew*exkurt*sigma^8*z^4 - 
             15552*beta*dz*skew*exkurt*sigma^8*z^2 - 6048*beta*dz*skew*exkurt*sigma^8 - 
             3456*beta*dz*skew*sigma^8*z^6 + 41472*beta*dz*skew*sigma^8*z^4 - 
             51840*beta*dz*skew*sigma^8*z^2 - 41472*beta*dz*skew*sigma^8 - 
             6912*beta*dz*exkurt*sigma^8*z^3 - 20736*beta*dz*sigma^8*z^3 - 
             1728*mu5*beta*dz*sigma^3*z^4 + 3456*mu5*beta*dz*sigma^3*z^2 + 
             1728*mu5*beta*dz*sigma^3 + 25*dz^2*skew^4*exkurt*sigma^8*z^12 - 
             450*dz^2*skew^4*exkurt*sigma^8*z^10 + 2475*dz^2*skew^4*exkurt*sigma^8*z^8 - 
             3900*dz^2*skew^4*exkurt*sigma^8*z^6 + 675*dz^2*skew^4*exkurt*sigma^8*z^4 + 
             1350*dz^2*skew^4*exkurt*sigma^8*z^2 + 225*dz^2*skew^4*exkurt*sigma^8 + 
             194*dz^2*skew^4*sigma^8*z^12 - 3252*dz^2*skew^4*sigma^8*z^10 + 
             16566*dz^2*skew^4*sigma^8*z^8 - 24024*dz^2*skew^4*sigma^8*z^6 + 
             3798*dz^2*skew^4*sigma^8*z^4 + 6876*dz^2*skew^4*sigma^8*z^2 + 
             1026*dz^2*skew^4*sigma^8 + 240*dz^2*skew^3*exkurt*sigma^8*z^9 - 
             2160*dz^2*skew^3*exkurt*sigma^8*z^7 + 2160*dz^2*skew^3*exkurt*sigma^8*z^5 + 
             720*dz^2*skew^3*exkurt*sigma^8*z^3 + 2016*dz^2*skew^3*sigma^8*z^9 - 
             16992*dz^2*skew^3*sigma^8*z^7 + 15840*dz^2*skew^3*sigma^8*z^5 + 
             4896*dz^2*skew^3*sigma^8*z^3 - 40*mu5*dz^2*skew^3*sigma^3*z^12 + 
             720*mu5*dz^2*skew^3*sigma^3*z^10 - 3960*mu5*dz^2*skew^3*sigma^3*z^8 + 
             6240*mu5*dz^2*skew^3*sigma^3*z^6 - 1080*mu5*dz^2*skew^3*sigma^3*z^4 - 
             2160*mu5*dz^2*skew^3*sigma^3*z^2 - 360*mu5*dz^2*skew^3*sigma^3 + 
             90*dz^2*skew^2*exkurt^2*sigma^8*z^10 - 990*dz^2*skew^2*exkurt^2*sigma^8*z^8 + 
             2340*dz^2*skew^2*exkurt^2*sigma^8*z^6 - 540*dz^2*skew^2*exkurt^2*sigma^8*z^4 - 
             1350*dz^2*skew^2*exkurt^2*sigma^8*z^2 - 270*dz^2*skew^2*exkurt^2*sigma^8 - 
             96*dz^2*skew^2*exkurt*sigma^8*z^12 + 2376*dz^2*skew^2*exkurt*sigma^8*z^10 - 
             16200*dz^2*skew^2*exkurt*sigma^8*z^8 + 29952*dz^2*skew^2*exkurt*sigma^8*z^6 + 
             864*dz^2*skew^2*exkurt*sigma^8*z^4 - 19656*dz^2*skew^2*exkurt*sigma^8*z^2 - 
             4536*dz^2*skew^2*exkurt*sigma^8 - 144*dz^2*skew^2*sigma^8*z^12 + 
             4500*dz^2*skew^2*sigma^8*z^10 - 32940*dz^2*skew^2*sigma^8*z^8 + 
             64296*dz^2*skew^2*sigma^8*z^6 + 19512*dz^2*skew^2*sigma^8*z^4 - 
             53964*dz^2*skew^2*sigma^8*z^2 - 12492*dz^2*skew^2*sigma^8 - 
             432*mu5*dz^2*skew^2*sigma^3*z^9 + 3888*mu5*dz^2*skew^2*sigma^3*z^7 - 
             3888*mu5*dz^2*skew^2*sigma^3*z^5 - 1296*mu5*dz^2*skew^2*sigma^3*z^3 + 
             16*mu6*dz^2*skew^2*sigma^2*z^12 - 348*mu6*dz^2*skew^2*sigma^2*z^10 + 
             2244*mu6*dz^2*skew^2*sigma^2*z^8 - 4056*mu6*dz^2*skew^2*sigma^2*z^6 + 
             792*mu6*dz^2*skew^2*sigma^2*z^4 + 1764*mu6*dz^2*skew^2*sigma^2*z^2 + 
             324*mu6*dz^2*skew^2*sigma^2 + 432*dz^2*skew*exkurt^2*sigma^8*z^7 - 
             864*dz^2*skew*exkurt^2*sigma^8*z^5 - 432*dz^2*skew*exkurt^2*sigma^8*z^3 - 
             1152*dz^2*skew*exkurt*sigma^8*z^9 + 13536*dz^2*skew*exkurt*sigma^8*z^7 - 
             16704*dz^2*skew*exkurt*sigma^8*z^5 - 10080*dz^2*skew*exkurt*sigma^8*z^3 - 
             72*mu5*dz^2*skew*exkurt*sigma^3*z^10 + 792*mu5*dz^2*skew*exkurt*sigma^3*z^8 - 
             1872*mu5*dz^2*skew*exkurt*sigma^3*z^6 + 432*mu5*dz^2*skew*exkurt*sigma^3*z^4 + 
             1080*mu5*dz^2*skew*exkurt*sigma^3*z^2 + 216*mu5*dz^2*skew*exkurt*sigma^3 - 
             1728*dz^2*skew*sigma^8*z^9 + 25920*dz^2*skew*sigma^8*z^7 - 
             36288*dz^2*skew*sigma^8*z^5 - 36288*dz^2*skew*sigma^8*z^3 - 
             432*mu5*dz^2*skew*sigma^3*z^10 + 4464*mu5*dz^2*skew*sigma^3*z^8 - 
             10656*mu5*dz^2*skew*sigma^3*z^6 - 3168*mu5*dz^2*skew*sigma^3*z^4 + 
             10512*mu5*dz^2*skew*sigma^3*z^2 + 2736*mu5*dz^2*skew*sigma^3 + 
             192*mu6*dz^2*skew*sigma^2*z^9 - 2016*mu6*dz^2*skew*sigma^2*z^7 + 
             2304*mu6*dz^2*skew*sigma^2*z^5 + 864*mu6*dz^2*skew*sigma^2*z^3 + 
             48*mu7*dz^2*skew*sigma*z^10 - 528*mu7*dz^2*skew*sigma*z^8 + 
             1248*mu7*dz^2*skew*sigma*z^6 - 288*mu7*dz^2*skew*sigma*z^4 - 
             720*mu7*dz^2*skew*sigma*z^2 - 144*mu7*dz^2*skew*sigma + 
             81*dz^2*exkurt^3*sigma^8*z^8 - 324*dz^2*exkurt^3*sigma^8*z^6 + 
             162*dz^2*exkurt^3*sigma^8*z^4 + 324*dz^2*exkurt^3*sigma^8*z^2 + 
             81*dz^2*exkurt^3*sigma^8 + 882*dz^2*exkurt^2*sigma^8*z^8 - 
             3528*dz^2*exkurt^2*sigma^8*z^6 + 468*dz^2*exkurt^2*sigma^8*z^4 + 
             6120*dz^2*exkurt^2*sigma^8*z^2 + 2178*dz^2*exkurt^2*sigma^8 + 
             3132*dz^2*exkurt*sigma^8*z^8 - 15984*dz^2*exkurt*sigma^8*z^6 - 
             2376*dz^2*exkurt*sigma^8*z^4 + 29808*dz^2*exkurt*sigma^8*z^2 + 
             16956*dz^2*exkurt*sigma^8 - 432*mu5*dz^2*exkurt*sigma^3*z^7 + 
             864*mu5*dz^2*exkurt*sigma^3*z^5 + 432*mu5*dz^2*exkurt*sigma^3*z^3 - 
             108*mu6*dz^2*exkurt*sigma^2*z^8 + 432*mu6*dz^2*exkurt*sigma^2*z^6 - 
             216*mu6*dz^2*exkurt*sigma^2*z^4 - 432*mu6*dz^2*exkurt*sigma^2*z^2 - 
             108*mu6*dz^2*exkurt*sigma^2 + 3564*dz^2*sigma^8*z^8 - 19440*dz^2*sigma^8*z^6 - 
             5832*dz^2*sigma^8*z^4 + 40176*dz^2*sigma^8*z^2 + 26892*dz^2*sigma^8 - 
             2592*mu5*dz^2*sigma^3*z^7 + 5184*mu5*dz^2*sigma^3*z^5 + 
             6048*mu5*dz^2*sigma^3*z^3 - 432*mu6*dz^2*sigma^2*z^8 + 
             2304*mu6*dz^2*sigma^2*z^6 - 3456*mu6*dz^2*sigma^2*z^2 - 
             1296*mu6*dz^2*sigma^2 + 288*mu7*dz^2*sigma*z^7 - 
             576*mu7*dz^2*sigma*z^5 - 288*mu7*dz^2*sigma*z^3 + 36*mu8*dz^2*z^8 - 
             144*mu8*dz^2*z^6 + 72*mu8*dz^2*z^4 + 144*mu8*dz^2*z^2 + 
             36*mu8*dz^2)/(20736*beta^2*sigma^6)
    
        names(est) <- names(variance) <- NULL
        list(est = est, se = sqrt(variance))
}