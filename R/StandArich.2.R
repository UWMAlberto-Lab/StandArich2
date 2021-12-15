# New function to run first the permutation  to estimate for all pops allelic richness and its standard deviation across loci
# for all samples sizes from n=1 to n=N (N=sample size) and second  tables with average allelic richness and its standard
# deviation across permutations of similar n

rgenotypes.arich<-function (Popdata, n.replicates,loc.labels){
  Tdata <- Popdata
  popnames<-unique(Tdata[,1])

  npops<-length(popnames)
  n.all <- length(Tdata[1, ]) - 2
  n.loc <- n.all/2
  i <- 1
  d.lim <- 1:npops

  repeat {
    c1 <- Tdata[, 1] == popnames[i]
    d.lim[i] <- length(Tdata[c1, 1])
    i <- i + 1
    if (i > npops)
      break
  }

  res.subsample <-data.frame(matrix(nrow=(n.replicates * sum(d.lim)),ncol=(n.loc +4)))
  names(res.subsample)<- c("Pop", "Ind", "A", "sdA", loc.labels)


  PopSize<-tapply(Tdata[,1],factor(Tdata[,1],levels=unique(Tdata[,1])),length)

  l<-1
  for(p in 1:npops){

    TdataPop<-Tdata[Tdata[,1]==popnames[p],]

    for(n in 1:PopSize[p]){
      for(r in 1:n.replicates){
        c.a1<-3
        c.a2<-4
        for(loc in 1:n.loc){
          repeat{
            indexInd<- sample(1:PopSize[p],n)
            AlleleBag<-c(TdataPop[indexInd,c.a1],TdataPop[indexInd,c.a2])
            AlleleBag<-AlleleBag[AlleleBag!=999]
            if(length(AlleleBag)>0)break
          }
          res.subsample[l,1]<-as.character(popnames[p])
          res.subsample[l,2]<-n
          res.subsample[l,(4+loc)] <-length(unique(AlleleBag))
          c.a1<-c.a1+2
          c.a2<-c.a2+2
        }
        l<-l+1
      }
    }

  }

  res.subsample$A<-apply(res.subsample[,-(1:4)],1,mean)
  res.subsample$sdA<-apply(res.subsample[,-(1:4)],1,sd)

  # Summarize the results in res.subsample. This replaces the need to use a separate function (ie standArich)
  StandArich.OUT<-list(
    R=res.subsample ,
    A=tapply(res.subsample$A,list(factor(res.subsample$Pop,levels=unique(res.subsample$Pop)),res.subsample$Ind),mean),
    A.sd=tapply(res.subsample$A,list(factor(res.subsample$Pop,levels=unique(res.subsample$Pop)),res.subsample$Ind),sd)
  )

  StandArich.OUT
}

##################################################### NEW allele.freq.plot###################################


plot.allele.freq<-function( data, LocusNames,Plogbase=20,textlegend=0.8,freqlegend=c(0.1,0.15,0.2,0.5,0.75,1),
                            psize=7,poptext=0.5,allelesize=1,xaxispos=0,yplotlim=-3,alleleangle=45,
                            pdfout=FALSE,pdfname="out.pdf"){


  PopCodes<-unique(data[,1])
  npop<-length(PopCodes)
  nloc<-(ncol(data)-2)/2
  a<-3

  if(pdfout==T | pdfout==TRUE){
    pdf(file=pdfname)}

  repeat{
    a1stcol<-data[,a]
    a2ndcol<-data[,a+1]
    AllelesL<-unique(c(a1stcol,a2ndcol))
    AllelesL<-AllelesL[AllelesL!=999]
    AllelesL<-sort(AllelesL)

    nallL<-length(AllelesL)

    dataL<-data[,c(1,a:(a+1))] #filtering the data by locus

    #starting the plot Sys.infor used to detect the operating system
    if(pdfout==F | pdfout==FALSE){
      switch(Sys.info()[['sysname']],
             Windows= {windows()},
             Linux  = {x11()},
             Darwin = {x11()})
    }

    plot(1,1,type="n",xlab="alleles",ylab="",axes=F,
         xlim=c(0,nallL+3),ylim=c(yplotlim,npop),main=LocusNames[(((a+1)/2)-1)],)
    axis(1,pos=xaxispos, at=seq(1,nallL),labels=rep("",nallL))
    text(seq(1,nallL),rep(-2,nallL),labels=AllelesL,cex=allelesize,srt=alleleangle)

    text(rep(0,npop),seq(npop,1,-1),PopCodes,cex=poptext)

    p<-1 #pop counter
    repeat{
      dataLPop<-dataL[dataL[,1]==PopCodes[p],] # Filtering data by pop, previously filtered by locus

      PopAll<-sort(unique(c(dataLPop[,2],dataLPop[,3])))
      PopAll<-PopAll[PopAll!=999]
      nPopAll<-length(PopAll)
      LPopAllCount<-tapply(c(dataLPop[,2],dataLPop[,3]),as.factor(c(dataLPop[,2],dataLPop[,3])),length)
      LPopAllCount<-LPopAllCount[names(LPopAllCount)!="999"]
      LPopAllFreq<-LPopAllCount/sum(LPopAllCount)

      for(ap in 1:nPopAll ){
        points(grep(PopAll[ap],AllelesL),(npop+1)-p,
               cex=(log(LPopAllFreq[ap]+1,base=Plogbase)*psize),pch=16)
      }

      p<-p+1
      if(p>npop)break}

    #legend
    text(rep(nallL+2,length(freqlegend)),seq(npop,npop-((length(freqlegend)-1)*2),-2),
         adj=0,cex=textlegend,labels=freqlegend)

    points(rep(nallL+1,length(freqlegend)),seq(npop,npop-((length(freqlegend)-1)*2),-2)
           ,cex=(log(freqlegend+1,base=Plogbase)*psize),pch=16)

    a<-a+2

    if(a>ncol(data))
      break}

  if(pdfout==T | pdfout==TRUE){
    dev.off()}
}


################################## Allele.genotype.plot
allele.genotype.plot<-function (results, g=0,xmin=0,xmax=50,xmark=10,main="",lwd=1,lty=1,lcol="black",print.pop=FALSE,pop.text.size=0.5)
{
    n.loci <- (length(results[1, ]) - 4)
    #attach(results)
    Population <- factor(results$Pop,levels=unique(results$Pop))
    Pops <- levels(Population)
    npops<-length(Pops)
    lty<-if(length(lty)>1){lty}else{rep(lty,npops)}
    if(length(lty)!=npops){stop("lty vector different from number of pops")}
    lwd<-if(length(lwd)>1){lwd}else{rep(lwd,npops)}
    if(length(lwd)!=npops){stop("lwd vector different from number of pops")}
    lcol<-if(length(lcol)>1){lcol}else{rep(lcol,npops)}
    if(length(lcol)!=npops){stop("lcol vector different from number of pops")}
        Sample.size <- results$Ind
    Mean.locus <- results$A
    var.locus <- results$sdA
    d.lim <- tapply(results[, 2], Population, max)
    max.y <- round(max(Mean.locus) + 1,digits=0)
    max.x <- round(max(as.numeric(Sample.size)),digits=0)

    min.size <- g
    i <- 1
    options(warn=-1)
#1:max.x, 1:max.y
    plot(1:max.x, type = "n", xlim = range(1, max.x + 10), ylim = range(1,
        max.y), axes = FALSE, xlab = "Nº of genotypes", ylab = "Allelic richness",main=main)
    axis(1, pos = 0.8, at = seq(xmin, xmax, xmark))
    axis(2, pos = 0, at = seq(1, max.y, 1), las = 1)
    repeat {
        ssize <- 1
        mean.ssize <- 1:d.lim[i]
        C1 <- Population == Pops[i]
        repeat {
            C2 <- Sample.size == ssize
            C3 <- C1 & C2
            mean.ssize[ssize] <- mean(Mean.locus[C3])
            ssize <- ssize + 1
            if (ssize > d.lim[i])
                break
        }
        lines(x=(1:d.lim[i]), y=mean.ssize,lwd=lcol[i],lty=lty[i])
        if(print.pop){text(x=(d.lim[i]+3),y=mean.ssize[d.lim[i]],Pops[i],cex=pop.text.size)}
        i <- i + 1
        if (i > length(d.lim))
            break
    }
   options(warn=0)
    lines(c(g,g),c(1,max.y),lty = 4)

}
