library(ggplot2)
library(data.table)
setwd("~/Git/k-similar-neighbor/")

plot_corr_matrix <- function(assoc_file, sub_ticks=NULL)
{
  
  X <- read.table(assoc_file, header=TRUE, row.names=1)
  tmp <- fread("data/SandP_companies.csv")
  tmp <- tmp[order(sector, subindustry, ticker)][ticker %in% names(X)]
  tmp[, idx := 1:.N]
  X <- X[tmp$ticker, tmp$ticker]
  X$ticks <- names(X)
  X$ticks <- factor(X$ticks, levels=X$ticks)
  X <- as.data.table(X)
  idxsec <- which(tmp[1:(.N - 1), sector] != tmp[2:.N, sector])
  idxlab <- tmp[, unique(sector)]
  idxloc <- tmp[, round(mean(idx)), by=sector][, V1]
  
  X <- melt(X, id.vars="ticks", measure.vars=setdiff(names(X), "ticks"))
  ggplot(data=X, aes(x=ticks, y=variable, fill=value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    #theme(axis.text.x=element_text(angle = -15)) +
    theme(axis.text=element_blank()) +
    annotate(geom="text", x=idxloc, y=idxloc, label=idxlab) +
    scale_x_discrete(breaks=tmp$ticker[idxloc], labels=idxlab) + 
    scale_y_discrete(breaks=tmp$ticker[idxloc], labels=idxlab) + 
    geom_vline(xintercept=idxsec, size=0.5, lty=3) + 
    geom_hline(yintercept=idxsec, size=0.5, lty=3)
}

plot_corr_matrix("data/association_results/corr_20190201_20200214_483")
