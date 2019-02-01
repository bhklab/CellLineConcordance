plotSubHeatmap <- function(x.m, levels=1){
  ## Just the match plot
  if(levels==1){
    p <- ggplot(x.m, aes(X2, X1)) + 
      geom_tile(aes(fill = value),colour = "black") +
      scale_fill_gradientn(limits=c(0.5, 1), breaks=seq(0.6,1,by=0.05), colours=sunset.col) +
      #scale_fill_gradient2(limits=c(0.7, 1), breaks=seq(0.7,1,by=0.05), low="black", mid="red", high="white", midpoint=0.85) +
      coord_fixed(ratio=1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
            axis.text.y = element_text(size=8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      facet_grid(fingerprint ~ anno)
    grid.arrange(p, ncol=2, nrow=2,
                 heights=c(2.5,2.5) )
  }   
  
  ## Both Match and Mismatch Plots
  if(levels==2){
    # Split to match and mismatch based on "anno" col
    split.melt <- split(x.m, x.m[,'anno'])
    
    # Match plot
    match.p <- ggplot(split.melt[['Matching annotation']], aes(X2, X1)) + 
      geom_tile(aes(fill = value),colour = "black") +
      scale_fill_gradientn(limits=c(0.5, 1), breaks=seq(0.6,1,by=0.05), colours=sunset.col, guide=FALSE) +
      #scale_fill_gradient2(limits=c(0.7, 1), breaks=seq(0.7,1,by=0.05), low="black", mid="red", high="white", midpoint=0.85) +
      coord_fixed(ratio = 1) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
            axis.text.y = element_text(size=8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(1, 0, 1, 1), "lines")) +
      facet_grid(. ~ anno)
    
    #Mismatch plot
    mismatch.p <- ggplot(split.melt[['Different annotation']], aes(X2, X1)) + 
      geom_tile(aes(fill = value),colour = "black") +
      scale_fill_gradientn(limits=c(0.5, 1), breaks=seq(0.6,1,by=0.05), colours=sunset.col) +
      #scale_fill_gradient2(limits=c(0.7, 1), breaks=seq(0.7,1,by=0.05), low="black", mid="red", high="white", midpoint=0.85) +
      coord_fixed(ratio = 1) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(1, 1, 1, 0), "lines")) +
      facet_grid(fingerprint ~ anno)
    
    grid.arrange(match.p, mismatch.p, ncol=2, nrow=2,
                 heights=c(2.5,2.5) )
  }
}
