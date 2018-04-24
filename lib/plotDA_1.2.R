plotDA <- function(scores, loadings=NULL, fac, groupfac=NULL, npoints=100, level=0.95,
                   facname="Factor", groupfacname=NULL, theme="theme_bw", coord.equal=FALSE,
                   propLoadings = 1,
                   xlab="First discriminant axis", ylab="Second discriminant axis") {
  
  # librairies
  require(ellipse)
  require(ggplot2)
  require(grid)
  require(plyr)
  
  colnames(scores) <- c("DS1", "DS2")
  if (!is.null(loadings)) colnames(loadings) <- c("DS1", "DS2")
  
  alpha <- 1-level

  centersFac <- data.frame(DS1=rep(NA, length(levels(fac))),
                           DS2=rep(NA, length(levels(fac))))
  rownames(centersFac) <- levels(fac)
  
  for (i in 1:length(levels(fac))) {
    centersFac[i,]<-apply(X=scores[fac == levels(fac)[i],], MARGIN=2, FUN="mean")
  }
  
  # error ellipses (conf. region about mean)
  ellErr<-list() 
  for (i in 1:length(levels(fac))) {
    tmp<-scores[fac == levels(fac)[i],]
    ellErr[[i]]<-ellipse(cor(tmp),
                         scale=c(sd(tmp[,1])/sqrt(nrow(tmp)), sd(tmp[,2])/sqrt(nrow(tmp))),
                         centre=apply(tmp, 2, "mean"),
                         level=level, npoints=npoints
    )
  }
  names(ellErr) <- levels(fac)
  
  # deviation ellipses (conf. region about population)
  ellDev<-list() 
  for (i in 1:length(levels(fac))) {
    tmp<-scores[fac == levels(fac)[i],]
    ellDev[[i]]<-ellipse(cor(tmp),
                         scale=c(sd(tmp[,1]), sd(tmp[,2])),
                         centre=c(mean(tmp[,1]),mean(tmp[,2])),
                         level=level, npoints=npoints
    )
  }
  names(ellDev) <- levels(fac)
  
  if (!is.null(groupfac)) {
    ellDevGrp<-list() # deviation ellipses (conf. region about population)
    for (i in 1:length(levels(groupfac))) {
      tmp<-scores[groupfac == levels(groupfac)[i],]
      ellDevGrp[[i]]<-ellipse(cor(tmp),
                              scale=c(sd(tmp[,1]), sd(tmp[,2])),
                              centre=c(mean(tmp[,1]),mean(tmp[,2])),
                              level=level, npoints=npoints
      )
    }
    names(ellDevGrp) <- levels(groupfac)
  }
  
  # ellipse peut entrer en conflit avec le package "car"
  detach(package:ellipse)
  
  # Ellipses dans des data.frames (grâce au package "plyr")
  dataEllErr <- ldply(.data = ellErr, data.frame)
  dataEllDev <- ldply(.data = ellDev, data.frame)
  if (!is.null(groupfac)) {
    dataEllDevGrp <- ldply(.data = ellDevGrp, data.frame)
  }

  #PLOT
  ggScores <- data.frame(fac, scores)
  ggCentersFac <- data.frame(fac=levels(fac), centersFac)
  ggMontage <- ggplot(data=ggScores, aes(x=DS1, y=DS2)) +
    geom_vline(xintercept = 0, colour='grey15', lwd=0.3) +
    geom_hline(yintercept = 0, colour='grey15', lwd=0.3) +
    geom_polygon(data=dataEllDev, aes(group=.id), fill="grey50", colour="grey70",
                 alpha=0.2)
  
  if(theme == "theme_bw") {
    ggMontage <- ggMontage + geom_point(aes(shape = fac), size=2, alpha=0.7)
  } else {
    ggMontage <- ggMontage + geom_point(aes(colour = fac), size=3, alpha=0.8)
  }
  
  ggMontage <- ggMontage + 
    geom_polygon(data=dataEllErr, aes(group=.id), fill="white",
                 alpha=0.6, colour="black") + #, alpha=0.8
    theme(axis.text.x=element_text(colour='black', size=11),
          axis.text.y=element_text(colour='black', size=11)) +
    xlab(xlab) + ylab(ylab) +
    geom_text(data=ggCentersFac, aes(x=DS1, y=DS2, label=fac), size=4,angle=45) +
    labs(colour = facname, shape=facname)
    
    
  if(!is.null(groupfac)) {
    ggMontage <- ggMontage + geom_polygon(data=dataEllDevGrp, aes(group=.id), fill=NA, alpha=0.8, colour='black')
  }
  
  if(!is.null(loadings)) {
    ggLoadings <- data.frame(ilrDef=rownames(loadings), propLoadings*loadings)
    ggMontage <- ggMontage + 
      geom_segment(data=ggLoadings, aes(x=0, y=0, xend=DS1, yend=DS2), size=0.5,
                   arrow=arrow(length=unit(0.2,"cm"))) +
      geom_text(data=ggLoadings, aes(DS1, DS2, label = ilrDef), hjust=0.5, vjust=-0.5, size=4)
  }
  
  if(coord.equal) {
    ggMontage <- ggMontage + coord_equal()
  }

  if(theme == "theme_bw") {
    ggMontage <- ggMontage + theme_bw()
  }
  
  ggMontage
}









