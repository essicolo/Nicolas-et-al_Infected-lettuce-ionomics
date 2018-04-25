plotDA <- function(discrimin, fac, groupfac=NULL, npoints=100, level=0.95,
                   scaling = 1,
                   circle_of_equilibrium = FALSE,
                   facname="Factor", 
                   groupfacname=NULL, 
                   propLoadings = 1,
                   xlab="First discriminant axis", 
                   ylab="Second discriminant axis") {
  
  # librairies
  library("ellipse")
  library("ggplot2")
  library("ggrepel")
  library("grid")
  library("plyr")
  
  if (scaling == 1) {
    scores  = discrimin$li
    loadings = discrimin$fa
  } else if (scaling == 2) {
    scores  = as.matrix(discrimin$li) %*% diag(discrimin$eig[1:2]**-0.5)
    loadings = as.matrix(discrimin$fa) %*% diag(discrimin$eig[1:2]**0.5)
  } else {
    print("scaling should be 1 or 2")
  }
  
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
    geom_polygon(data=dataEllDev, aes(group=.id), fill="grey50", colour="grey70", alpha=0.2) +
    geom_point(aes(shape = fac), alpha=0.6, fill='white', size=3) +
    geom_polygon(data=dataEllErr, aes(group=.id), fill="white",
                 alpha=0.6, colour="black") + #, alpha=0.8
    theme(axis.text.x=element_text(colour='black', size=11),
          axis.text.y=element_text(colour='black', size=11)) +
    xlab(xlab) + ylab(ylab) +
    geom_label_repel(data=ggCentersFac, aes(x=DS1, y=DS2, label=fac), size=4,
               fill='black', colour='white') +
    labs(fill = facname, shape=facname)
  
  
  if(!is.null(groupfac)) {
    ggMontage <- ggMontage + geom_polygon(data=dataEllDevGrp, aes(group=.id), fill=NA, alpha=0.8, colour='black')
  }
  
  if(circle_of_equilibrium & scaling == 1) {
    r = sqrt(2/nrow(loadings)) * propLoadings
    tt <- seq(0,2*pi,length.out = npoints)
    cex <- r * cos(tt)
    cey <-  r * sin(tt)
    ggMontage <- ggMontage + geom_path(data=data.frame(cex, cey), aes(cex, cey))
  }
  
  ggLoadings <- data.frame(ilrDef=rownames(loadings), 
                           propLoadings*loadings,
                           hjust = ifelse(loadings[, 1] < 0, 1.05, -0.05),
                           vjust = ifelse(loadings[, 2] < 0, 1.05, -0.05))
  
  ggMontage <- ggMontage + 
    geom_segment(data=ggLoadings, aes(x=0, y=0, xend=DS1, yend=DS2), size=0.5,
                 arrow=arrow(length=unit(0.2,"cm"))) +
    #geom_label(data=ggLoadings, aes(DS1, DS2, label = ilrDef, hjust=hjust, vjust=vjust), 
    #           size=4, fill='black', colour='white')
    geom_label_repel(data=ggLoadings, aes(DS1, DS2, label = ilrDef, hjust=hjust, vjust=vjust), 
               size=4, segment.size=0.2)
  
  ggMontage
}









