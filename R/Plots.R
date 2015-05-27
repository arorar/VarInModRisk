colors <- c("red","blue", "green", "brown", "orange", "black", "grey", "purple")

plot.efficiency <- function(distribution, method, size = c(200, 500, 1000), 
                            etl=FALSE, plot = TRUE) {
    
    temp.func <- function(se=TRUE) {
        main.grid <- if (method == "large-sample" && distribution == "gaussian")
            risk.largesample.gaussian.efficiency(etl) 
        else if (method == "large-sample" && distribution == "t")
            risk.largesample.t.efficiency(etl)
        else if (method == "large-sample" && distribution == "skew-t")
            risk.largesample.skewt.efficiency(etl)
        else if (method == "small-sample" && distribution == "gaussian")
            risk.smallsample.gaussian.efficiency(etl, size = size)
        else if (method == "small-sample" && distribution == "t")
            risk.smallsample.t.efficiency(etl, size = size)
        else if (method == "small-sample" && distribution == "skew-t")
            risk.smallsample.skewt.efficiency(etl, size = size)
        
        
        drop.param <- if(se) "est.bias" else "se.efficiency"            
        drop.param <- which(colnames(main.grid) == drop.param)
        data <- main.grid[,-drop.param]
        colnames(data) <-   
            if (method == "large-sample") {
                if (distribution == "gaussian") c("sign","param") 
                else if (distribution == "t") c("df","sign","param")
                else if (distribution == "skew-t") c("slant","df","sign","param")
            } else if (method == "small-sample") {
                if (distribution == "gaussian") c("size","sign","param") 
                else if (distribution == "t") c("size", "df", "sign", "param") 
                else if (distribution == "skew-t") c("size","slant","df","sign","param")
            }
        
        data <- na.omit(data)
        
        y.lab <- if (se && !etl) expression(paste("%",frac(se(VaR),se(mVaR))))
        else if (!se && !etl) "% Relative Bias"
        else if (se &&   etl) expression(paste("%",frac(se(ES),se(mES))))
        else if (!se &&  etl) "% Relative Bias"
        
        
        fivep.label <- "5%"; onep.label <- "1%"        
        
        p <- if (method == "large-sample") {
            if (distribution == "gaussian")
                ggplot(data, aes(x=sign,y=param)) + geom_line() + ylab(y.lab) + 
                xlab(bquote(alpha))
            else if (distribution == "t") {
                sign.filter <- c(0.01,0.05)
                data <- data[data[,"sign"] %in% sign.filter,]
                ggplot(data, aes(x=df,y=param)) + 
                    geom_line(aes(color=as.factor(sign))) +
                    scale_colour_manual(values = c("0.01"="red","0.05"="blue"), 
                            labels=c("1%", "5%"), name="Tail Probability") +
                    ylab(y.lab) + xlab(bquote(nu))
            } else if (distribution == "skew-t") {
                sign.filter <- c(0.01,0.05)
                data <- data[data[,"sign"] %in% sign.filter,]  
                
                ggplot(data, aes(x=df,y=param)) + facet_grid(sign~slant) + 
                    geom_line(aes(color=as.factor(sign))) +
                    scale_colour_manual(values = c("0.01"="red","0.05"="blue"), 
                                        labels=c("1%", "5%"), name="Tail Probability") +
                    scale_x_continuous(breaks=c(10,20,30)) + 
                    scale_y_continuous(breaks=pretty_breaks(n=4)) + 
                    ylab(y.lab) + xlab(bquote(nu))
                
            }
            
        } else if (method == "small-sample") {  
            
            size.col <- colors[1:length(size)]
            names(size.col) <- size
            
            if (distribution == "gaussian")
                ggplot(data, aes(x=sign,y=param)) + ylab(y.lab) + 
                xlab(bquote(alpha)) + geom_line(aes(color=as.factor(size))) + 
                scale_colour_manual(values = size.col, 
                                    labels = as.character(size), 
                                    name="Sample Size")
            else if (distribution == "t") {
                sign.filter <- c(0.01,0.05)
                data <- data[data[,"sign"] %in% sign.filter,]
                
                ggplot(data, aes(x=df,y=param)) + facet_grid(.~sign) + 
                    ylab(y.lab) + xlab(bquote(nu)) + 
                    geom_smooth(aes(color=as.factor(size)), se=FALSE, method=loess) + 
                    scale_colour_manual(values = size.col, labels=as.character(size), 
                                        name="Sample Size")
            } else if (distribution == "skew-t") {
                sign.filter <- c(0.01,0.05)
                data <- data[data[,"sign"] %in% sign.filter,]             
               
                ggplot(data, aes(x=df,y=param)) + facet_grid(sign~slant) + 
                    ylab(y.lab) + xlab(bquote(nu)) + 
                    geom_smooth(aes(color=as.factor(size)), se=FALSE, method=loess) + 
                    scale_x_continuous(breaks=c(10,20,30)) + 
                    scale_y_continuous(breaks=pretty_breaks(n=4)) + 
                    scale_colour_manual(values = size.col, labels=as.character(size), 
                                        name="Sample Size")
            }
        }
        
        p
    }
    
    risk <- if(etl) "mES" else "mVaR"
    title <- paste(risk, " bias and efficiency for ", 
                   ifelse(distribution == "gaussian", "Normal", 
                          ifelse(distribution == "t", "Standardized t","Standardized Skew-t")),
                   " distribution")
    p1 <- temp.func(se = FALSE); p2 <- temp.func()
    
    p <- if(distribution=="gaussian" && method == "large-sample") 
        arrangeGrob(p1, p2, nrow=1,ncol=2, main = textGrob(title,
                                       gp=gpar(fontface = "bold",fontsize=20)))
        else {
            g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs
            legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
            lheight <- sum(legend$height)
            arrangeGrob(
                arrangeGrob(p1 + theme(legend.position="none"), 
                            p2 + theme(legend.position="none"),
                            nrow=1),
                legend,
                ncol=1,
                heights = unit.c(unit(1, "npc") - lheight, lheight),
                main = textGrob(title, gp=gpar(fontface = "bold",fontsize=20)))
        }
            
    if (plot) print(p)
    p
}


