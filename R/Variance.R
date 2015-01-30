
factor <- function(z, distribution, param = NA, C = NA) {
    if(distribution == "gaussian") dnorm(z)
    else if(distribution == "t") {
        if ("nu" %in% names(param))  { 
            nu <- param$nu
            dt(sqrt(C)*z, df = nu - 2)/sqrt(C)
        } else NA
    }
}

# var.type should take care of which quantile to use
# distribution takes care of the moments
risk.ordinary <- function(alpha, method, distribution, param, quantile, moments, 
                          etl, C) {    

    z   <- quantile; dz <- factor(z, distribution, param, C = C)
    mu  <- moments["mu"];  mu2 <- moments["mu2"]
    mu3 <- moments["mu3"]; mu4 <- moments["mu4"]
    
    sigma   <- sqrt(mu2) 
    skew    <- mu3/mu2^(3/2) 
    exkurt  <- mu4/mu2^2 - 3
    
    sigma <- sqrt(mu2); skew <- mu3/mu2^(3/2); exkurt <- mu4/mu2^2 - 3
    
    est <- if(!etl) mu + C*sigma*z else mu - C*sigma/alpha*dz        
    
    variance <- 
        if(!etl) 
            ((C*z)^2/2*(1+exkurt/2) + skew*C*z + 1 )*sigma^2 
        else
            (((C*dz)^2*(exkurt + 2))/(4*alpha^2) - (C*dz*skew)/alpha + 1)*sigma^2 
    
    names(est) <- names(variance) <- NULL
    
    list(est = est, se = sqrt(variance))
}

# var.type should take care of which quantile to use
# distribution takes care of the moments
risk.modified <- function(alpha, method, distribution, quantile, moments, etl) {

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
        term1 <-  (dz*z^2)
        term2 <-  (z^4 - 2*z^2 - 1)*dz
        term3 <-  (z^6 - 9*z^4 + 9*z^2 + 3)*dz        
        
        factor <- dz + (term1)*skew/6 + (term2)*exkurt/24 + (term3)*skew^2/72
        metl <- mu - sigma/alpha*factor
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
        (sigma^8*(dz^2*(225*skew^4*exkurt + 1026*skew^4 - 270*skew^2*exkurt^2 - 
        4536*skew^2*exkurt - 12492*skew^2 + 81*exkurt^3 + 2178*exkurt^2 + 
        16956*exkurt + 26892) + dz^2*z^4*(675*skew^4*exkurt + 3798*skew^4 - 
        540*skew^2*exkurt^2 + 1440*skew^2*exkurt + 24696*skew^2 + 162*exkurt^3 + 
        468*exkurt^2 - 5832*exkurt - 11016) + dz^2*z^8*(2475*skew^4*exkurt + 
        16566*skew^4 - 990*skew^2*exkurt^2 - 16200*skew^2*exkurt - 32940*skew^2 + 
        81*exkurt^3 + 882*exkurt^2 + 3132*exkurt + 3564) - 
        dz^2*z^6*(3900*skew^4*exkurt + 24024*skew^4 - 2340*skew^2*exkurt^2 - 
        29376*skew^2*exkurt - 59112*skew^2 + 324*exkurt^3 + 3528*exkurt^2 + 
        12528*exkurt + 14256) + dz^2*z^2*(1350*skew^4*exkurt + 6876*skew^4 - 
        1350*skew^2*exkurt^2 - 19656*skew^2*exkurt - 53964*skew^2 + 
        324*exkurt^3 + 6120*exkurt^2 + 29808*exkurt + 40176) - 
        dz^2*skew^2*z^12*(96*exkurt - 25*skew^2*exkurt - 194*skew^2 + 144) + 
        6*dz^2*skew^2*z^10*(15*exkurt^2 - 542*skew^2 - 75*skew^2*exkurt + 
        396*exkurt + 750)) + 36*dz^2*mu8 + 20736*alpha^2*sigma^8 - 
        dz^2*z^6*(144*mu8 + 2592*mu5*sigma^3 - 288*mu7*sigma + 432*exkurt*mu5*sigma^3 - 
        3888*skew^2*mu5*sigma^3) + dz^2*z^2*(144*mu8 + 6048*mu5*sigma^3 - 
        288*mu7*sigma + 432*exkurt*mu5*sigma^3 - 1296*skew^2*mu5*sigma^3) + 
        dz^2*z^4*(72*mu8 + 5184*mu5*sigma^3 - 576*mu7*sigma + 
        864*exkurt*mu5*sigma^3 - 3888*skew^2*mu5*sigma^3) - 
        288*alpha*dz*(72*sigma^8*z^2 - 6*mu5*sigma^3 - 
        15*skew^3*sigma^7*sigma + 24*exkurt*sigma^8*z^2 + 
        144*skew*sigma^7*sigma - 24*skew^2*sigma^8*z^2 - 12*mu5*z^2*sigma^3 + 
        6*mu5*z^4*sigma^3 - 45*skew^3*sigma^7*z^2*sigma + 
        45*skew^3*sigma^7*z^4*sigma - 5*skew^3*sigma^7*z^6*sigma + 
        21*skew*exkurt*sigma^7*sigma + 180*skew*sigma^7*z^2*sigma - 
        144*skew*sigma^7*z^4*sigma + 12*skew*sigma^7*z^6*sigma + 
        54*skew*exkurt*sigma^7*z^2*sigma - 45*skew*exkurt*sigma^7*z^4*sigma + 
        4*skew*exkurt*sigma^7*z^6*sigma) + dz^2*z^8*(36*mu8 - 
        432*skew^2*mu5*sigma^3) - 4*dz^2*mu6*sigma^2*(27*exkurt + 
        108*exkurt*z^2 + 54*exkurt*z^4 - 108*exkurt*z^6 + 27*exkurt*z^8 - 
        81*skew^2 + 864*z^2 - 144*z^4 - 432*z^6 + 108*z^8 - 441*skew^2*z^2 - 
        198*skew^2*z^4 + 1014*skew^2*z^6 - 561*skew^2*z^8 + 87*skew^2*z^10 -
        4*skew^2*z^12 + 324) - 48*dz^2*skew*sigma*(3*mu7 + 15*mu7*z^2 +
        6*mu7*z^4 - 26*mu7*z^6 + 11*mu7*z^8 - mu7*z^10 - 18*mu6*z^2*sigma - 
        48*mu6*z^4*sigma + 42*mu6*z^6*sigma - 4*mu6*z^8*sigma) - 
        8*dz^2*skew*mu5*sigma^3*(234*exkurt*z^6 - 135*exkurt*z^2 - 
        54*exkurt*z^4 - 27*exkurt - 99*exkurt*z^8 + 9*exkurt*z^10 + 
        45*skew^2 - 1314*z^2 + 540*z^4 + 1188*z^6 - 558*z^8 + 54*z^10 + 
        270*skew^2*z^2 + 135*skew^2*z^4 - 780*skew^2*z^6 + 495*skew^2*z^8 - 
        90*skew^2*z^10 + 5*skew^2*z^12 - 342) - 
        48*dz^2*skew*sigma^7*z^2*sigma*(45*skew^2*exkurt*z^4 - 
        5*skew^2*exkurt*z^6 - 45*skew^2*exkurt*z^2 - 15*skew^2*exkurt - 
        42*skew^2*z^6 + 354*skew^2*z^4 - 330*skew^2*z^2 - 102*skew^2 - 
        9*exkurt^2*z^4 + 18*exkurt^2*z^2 + 9*exkurt^2 + 24*exkurt*z^6 - 
        282*exkurt*z^4 + 348*exkurt*z^2 + 210*exkurt + 36*z^6 - 540*z^4 + 
        756*z^2 + 756))/(20736*alpha^2*sigma^6)
        
        names(est) <- names(variance) <- NULL
        list(est = est, se = sqrt(variance))
}