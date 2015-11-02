colors <- c("red","blue", "green", "brown", "orange", "black", "grey", "purple")

plot.efficiency <- function(distribution, method, size = c(200, 500, 1000), 
                            replicates = 2500,
                            etl=FALSE, plot = TRUE) {
    
    temp.func <- function(se=TRUE) {
        main.grid <- if (method == "large-sample" && distribution == "gaussian")
            risk.largesample.gaussian.efficiency(etl) 
        else if (method == "large-sample" && distribution == "t")
            risk.largesample.t.efficiency(etl)
        else if (method == "small-sample" && distribution == "gaussian")
            risk.smallsample.gaussian.efficiency(etl, size = size, replicates = replicates)
        else if (method == "small-sample" && distribution == "t")
            risk.smallsample.t.efficiency(etl, size = size, replicates = replicates)

        drop.param <- if(se) "est.bias" else "se.efficiency"            
        drop.param <- which(colnames(main.grid) == drop.param)
        data <- main.grid[,-drop.param]
        colnames(data) <-   
            if (method == "large-sample") {
                if (distribution == "gaussian") c("sign","param") 
                else if (distribution == "t") c("df","sign","param")
            } else if (method == "small-sample") {
                if (distribution == "gaussian") c("size","sign","param") 
                else if (distribution == "t") c("size", "df", "sign", "param") 
            }
        
        data <- na.omit(data)
        
        y.lab <- if (se && !etl) "% EFF-mVaR"
        else if (!se && !etl) "% Relative Bias"
        else if (se &&   etl) "% EFF-mES"
        else if (!se &&  etl) "% Relative Bias"
        
        
        p <- if (method == "large-sample") {
            if (distribution == "gaussian")
                ggplot(data, aes(x=sign,y=param)) + geom_line() + ylab(y.lab) + 
                xlab(bquote(gamma))
            else if (distribution == "t") {
                sign.filter <- c(0.01, 0.025, 0.05)
                data <- data[data[,"sign"] %in% sign.filter,]
                ggplot(data, aes(x=df,y=param)) + 
                    geom_line(aes(color=as.factor(sign))) +
                    scale_colour_manual(values = c("0.01"="red", "0.025"="blue", "0.05"="green"), 
                            labels=c("1%", "2.5%", "5%"), name="Tail Probability") +
                    ylab(y.lab) + xlab(bquote(nu))
            } 
        } else if (method == "small-sample") {  
            
            size.col <- colors[1:length(size)]
            names(size.col) <- size
            
            if (distribution == "gaussian")
                ggplot(data, aes(x=sign,y=param)) + ylab(y.lab) + 
                xlab(bquote(gamma)) + geom_line(aes(color=as.factor(size))) + 
                scale_colour_manual(values = size.col, 
                                    labels = as.character(size), 
                                    name="Sample Size")
            else if (distribution == "t") {
                sign.filter <- c(0.01, 0.025, 0.05)
                data <- data[data[,"sign"] %in% sign.filter,]
                
                ggplot(data, aes(x=df,y=param)) + facet_grid(.~sign) + 
                    ylab(y.lab) + xlab(bquote(nu)) + 
                    geom_smooth(aes(color=as.factor(size)), se=FALSE, method=loess) + 
                    scale_colour_manual(values = size.col, labels=as.character(size), 
                                        name="Sample Size")
            } 
        }
        
        p
    }
    
    risk <- if(etl) "mES" else "mVaR"
    title <- paste(risk, " Bias and Efficiency for ", 
                   ifelse(distribution == "gaussian", "Normal", 
                          ifelse(distribution == "t", "Student's t")),
                   " Distribution")
    p1 <- temp.func(se = FALSE); p2 <- temp.func()
    
    p <- if(distribution=="gaussian" && method == "large-sample") 
        arrangeGrob(p1, p2, nrow=1,ncol=2, top = textGrob(title,
                                       gp=gpar(fontface = "bold",fontsize=18)))
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
                heights = unit.c(unit(0.9, "npc") - lheight, lheight),
                top = textGrob(title, gp=gpar(fontface = "bold",fontsize=18)))
        }
            
    if (plot) plot(p)
    p
}


