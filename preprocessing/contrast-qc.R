library(ggplot2)
library(xtable)
dir <- '/Users/rquevedo/Desktop/bhk_lab/data/hapmap3_normalized/apt_geno_qc'
setwd(dir)
dir.files <- list.files(dir, pattern = "*results.txt")

dataset.list <- list("gdsc" = read.csv(dir.files[3], header=T, sep="\t", comment.char="#"),
                     "cgp" = read.csv(dir.files[2], header = T, sep = "\t", comment.char = "#"),
                     "ccle" = read.csv(dir.files[1], header = T, sep = "\t", comment.char = "#"),
                     "pfizer" = read.csv(dir.files[4], header = T, sep = "\t", comment.char = "#"),
                     "gray" = read.csv(dir.files[5], header = T, sep = "\t", comment.char = "#"))

dataset.df <- as.data.frame(matrix(ncol=0, nrow=5,
                                   dimnames=list(c("GDSC", "CGP", "CCLE", "Pfizer", "Gray"), c())))

### Generate a summary table in pdf and tex format
dataset.df[,'n'] <- sapply(dataset.list, function(x) length(x[,'contrast.qc']))
dataset.df[,'mean-cqc'] <- sapply(dataset.list, function(x) mean(x[,'contrast.qc']))
dataset.df[,'median-cqc'] <- sapply(dataset.list, function(x) median(x[,'contrast.qc']))
dataset.df[,'low cqc*'] <- sapply(dataset.list, function(x) length(x[(x[,'contrast.qc'] <= 0.4),'contrast.qc']))
dataset.df[,'RE cqc ratio**'] <- sapply(dataset.list, function(x){
                                        high.ratio.col <- abs(x[,'contrast.qc.nsp'] - x[,'contrast.qc.sty']) >= 2
                                        return(length(x[high.ratio.col, 'contrast.qc']))
                                        })
outdir <- '/Users/rquevedo/Desktop/bhk_lab/results/qc/'
setwd(outdir)

#https://cloudconvert.com/tex-to-pdf
italic <- function(x){
  paste0('{\\emph{',x,'}}')
}
print(xtable(dataset.df, include.rownames=TRUE, floating=FALSE),
             sanitize.rownames.function= italic,
             booktabs=TRUE,
             file="/Users/rquevedo/Desktop/bhk_lab/results/qc/cqc_table.tex")
pdf("/Users/rquevedo/Desktop/bhk_lab/results/qc/cqc_table.pdf")
grid.table(dataset.df)
dev.off()

# Label and combine all the dataframes together into one dataframe
dataset.list <- lapply(names(dataset.list), function(x){
  dataset.rep <- rep(x, dim(dataset.list[[x]])[1])
  x <- cbind(x, dataset.list[[x]])
})
all.dataset.df <- do.call("rbind", dataset.list)

# Create ggplots for CQC values
all.cqc.df <- all.dataset.df[,c("x", "cel_files", "contrast.qc")]
all.cqc.m <- melt(all.cqc.df)
colnames(all.cqc.m) <- c("dataset", "cel_files", "metric", "value")

pcqc <- ggplot(all.cqc.m, aes(x=value)) +
  geom_density() +
  facet_grid(dataset ~ .)
pcqc.f <- pcqc + theme_bw() +
  annotate("rect", xmin=-Inf, xmax=0.4, ymin=0, ymax=Inf, fill="black", alpha=0.3) +
  scale_fill_manual(values = alpha(c("black"), .1)) +
  labs(x="", y="") +
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_blank())


# Create ggplots for Nsp Sty Ratio
all.dataset.df$nsp_sty_ratio <- abs(all.dataset.df$contrast.qc.nsp - all.dataset.df$contrast.qc.sty)
all.ratio.df <- all.dataset.df[,c("x", "cel_files", "nsp_sty_ratio")]
all.ratio.m <- melt(all.ratio.df)
colnames(all.ratio.m) <- c("dataset", "cel_files", "metric", "value")

pratio <- ggplot(all.ratio.m, aes(x=value)) +
  geom_density() +
  facet_grid(dataset ~ .)
pratio.f <- pratio + theme_bw() +
  labs(x="", y="") +
  theme(legend.position = "none") +
  annotate("rect", xmin=2.0, xmax=Inf, ymin=0, ymax=Inf, fill="black", alpha=0.3)

pdf("cqc_plots.pdf")
multiplot(pcqc.f, pratio.f, cols=2)
dev.off()

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
  
