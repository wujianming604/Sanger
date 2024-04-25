#!/usr/bin/env Rscript
library(sangerseqR)
library(stringr)
args=commandArgs(T)

name <- args[1]
localPos <- args[2]
pngName <- args[3]

lis <- as.character(unlist(strsplit(name, split = "-")))
#pngName <- paste0(lis[2],"_",lis[3],".png")
seq <- readsangerseq(name)
upStream_length <- 20
downStream_length <- 20

if(str_detect(localPos,',')){
  print("截取多个位点")
  start_pos <- as.numeric(strsplit(localPos, split=',')[[1]][1])
  end_pos <- as.numeric(strsplit(localPos, split=',')[[1]][2])
  mul <- end_pos - start_pos
  setMethod("chromatogram", "sangerseq", 
            function(obj, trim5=0, trim3=0, 
                    showcalls=c("primary", "secondary", "both", "none"), 
                    width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
                    filename=NULL, showtrim=FALSE, showhets=TRUE) {
              
    if (start_pos + downStream_length > nchar(seq@primarySeq)){
      stop("截图长度不能超过原图总长度！")
    }
    if (end_pos < upStream_length){
      stop("上游序列长度不能超过截取位置！")
    }
    
    upStream_start <- start_pos - upStream_length
    upStream_end <- start_pos -1
    downStream_start <- end_pos + 1
    downStream_end <- end_pos + downStream_length
    
    
    upStream_seq <- substr(seq@primarySeq,upStream_start,upStream_end)
    downStream_seq <- substr(seq@primarySeq,downStream_start,downStream_end)
    
    pattern <- paste0(upStream_seq,'([ACTG]*)',downStream_seq)
    
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    standSeq <- paste(basecalls1,collapse="")
    match <- str_extract_all(standSeq,pattern)
  
    
    revCom_up <- DNAStringSet(upStream_seq,)%>% reverse %>% complement
    revCom_down <- DNAStringSet(downStream_seq,)%>% reverse %>% complement
    revCom_pattern <- paste0(revCom_down,'([ACTG]*)',revCom_up)
    revCom_match  <- str_extract_all(standSeq,revCom_pattern)
    
    #target_seq <- str_extract_all(basecalls1,)
    #pic_str <- "CCTCCGTAAATACTGGACCCAAGTTACTGC" # wait_for_change
    pic_pos <- nchar(upStream_seq)  # wait_for_change
    originalpar <- par(no.readonly=TRUE)
    #showcalls <- showcalls[3]
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    list_str_basecalls1 <- paste(basecalls1,collapse="")
    list_tmp_store_1 <- strsplit(list_str_basecalls1, match)
    list_tmp_store_2 <- strsplit(list_str_basecalls1, revCom_match)
    
    
    if (length(list_tmp_store_1[[1]]) == 2)
    {
      trim5 = nchar(list_tmp_store_1[[1]][1]) 
      trim3 = nchar(list_tmp_store_1[[1]][2])
    }
    else if (length(list_tmp_store_2[[1]]) == 2)
    {
      trim5 = nchar(list_tmp_store_2[[1]][1])
      trim3 = nchar(list_tmp_store_2[[1]][2])
      print(trim5)
      print(trim3)
      #pic_pos <- nchar(pic_str) - pic_pos + 1
    }
    aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(aveposition)] 
    basecalls2 <- basecalls2[1:length(aveposition)] 
    if(showtrim == FALSE) {
      if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
      else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
      if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
      else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
      aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
    }
    indexes <- 1:length(basecalls1)
    trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
    #false if not trimmed
    if (!is.null(trim3)) {
      traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                              nrow(traces))), ]
    }
    if (!is.null(trim5)) {
      offset <- max(c(1, aveposition[1] - 10))
      traces <- traces[offset:nrow(traces),]
      aveposition <- aveposition - (offset-1)
    }
    maxsignal <- apply(traces, 1, max)
    ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
    p <- c(0, aveposition, nrow(traces))
    midp <- diff(p)/2
    starts <- aveposition - midp[1:(length(midp)-1)]
    starthets <- starts
    #print(starthets)
    for(ii in 1:length(starthets))
    {
      if(mul == 1){
        if(ii != pic_pos+1 & ii != pic_pos +2)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 2){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 3){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 4){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 5){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 6){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 7){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 9){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 10){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10 & ii != pic_pos +11)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 11){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10 & ii != pic_pos +11  & ii != pic_pos +12)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 15){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10 & ii != pic_pos +11  & ii != pic_pos +12  & ii != pic_pos +13  & ii != pic_pos +14  & ii != pic_pos +15 & ii != pic_pos +16)
        {
          starthets[ii] <- NA
        }
      }

    }
    ends <- aveposition + midp[2:(length(midp))]
    endhets <- ends
    for(ii in 1:length(endhets))
    {
      if(mul == 1){
        if(ii != pic_pos+1 & ii != pic_pos +2)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 2){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 3){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 4){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 5){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 6){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 7){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 9){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 10){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10 & ii != pic_pos +11)
        {
          starthets[ii] <- NA
        }
      }
      else if(mul == 11){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10 & ii != pic_pos +11  & ii != pic_pos +12)
        {
          endhets[ii] <- NA
        }
      }
      else if(mul == 15){
        if(ii != pic_pos+1 & ii != pic_pos +2 & ii != pic_pos +3 & ii != pic_pos +4  & ii != pic_pos +5  & ii != pic_pos +6  & ii != pic_pos +7  & ii != pic_pos +8 & ii != pic_pos +9 & ii != pic_pos +10 & ii != pic_pos +11  & ii != pic_pos +12  & ii != pic_pos +13  & ii != pic_pos +14  & ii != pic_pos +15 & ii != pic_pos +16)
        {
          endhets[ii] <- NA
        }
      }
    }
    starttrims <- starts
    starttrims[!trimmed] <- NA
    endtrims <- ends
    endtrims[!trimmed] <- NA
    
    colortranslate <- c(A="green", C="blue", G="black", T="red")
    colorvector1 <- unname(colortranslate[basecalls1])
    colorvector1[is.na(colorvector1)] <- "purple"
    colorvector2 <- unname(colortranslate[basecalls2])
    colorvector2[is.na(colorvector2)] <- "purple"
    
    valuesperbase <- nrow(traces)/length(basecalls1)
    tracewidth <- width*valuesperbase
    breaks <- seq(1,nrow(traces), by=tracewidth)
    numplots <- length(breaks)
    if(!is.null(filename)) png(filename, width=2500, height=600, pointsize = 8, res=500) 
    #par(mar=c(2,2,2,1), mfrow=c(numplots, 1))
    par(mar=c(0,0,1,0), mfrow=c(numplots, 1))
    basecallwarning1 = 0
    basecallwarning2 = 0
    j = 1
    
    for(i in breaks) {
      range <- aveposition >= i & aveposition < (i+tracewidth)
      starthet <- starthets[range] - tracewidth*(j-1)
      starthet[starthet < 0] <- 0
      endhet <- endhets[range] - tracewidth*(j-1)
      endhet[endhet > tracewidth] <- tracewidth
      lab1 <- basecalls1[range]
      lab2 <- basecalls2[range]
      pos <- aveposition[range] - tracewidth*(j-1)
      colors1 <- colorvector1[range]
      colors2 <- colorvector2[range]
      starttrim <- starttrims[range] - tracewidth*(j-1)
      endtrim <- endtrims[range] - tracewidth*(j-1)
      plotrange <- i:min(i+tracewidth, nrow(traces))
      plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", bty="n", xlab="", yaxt="n", xlim=c(23,tracewidth-23))
      #title(main = "Stopping Distance versus Speed", font.main= 4) # wait_for_change
      if (showhets==TRUE) {
        rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')#lty=3
      }
      if (showtrim==TRUE) {
        rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent', 
            density=15)
      }
      lines(traces[plotrange,1], col="green")
      lines(traces[plotrange,2], col="blue")
      lines(traces[plotrange,3], col="black")
      lines(traces[plotrange,4], col="red")
      #mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
      
      for(k in 1:length(lab1)) {
        if (showcalls=="primary" | showcalls=="both") {
          if (is.na(basecalls1[1]) & basecallwarning1==0) {
            warning("Primary basecalls missing")
            basecallwarning1 = 1
          } 
          else if (length(lab1) > 0) {   
            axis(side=3, at=pos[k], labels=lab1[k], cex.axis=1, col.axis=colors1[k], 
                family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0, 
                                                          -1), tick=FALSE)
          }
        }
        if (showcalls=="secondary" | showcalls=="both") {
          if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
            warning("Secondary basecalls missing")
            basecallwarning2 = 1
          } 
          else if (length(lab2) > 0) { 
            axis(side=3, at=pos[k], labels=lab2[k], cex.axis=1, col.axis=colors2[k], 
                family="mono", cex=cex.base, line=-1.0, tick=FALSE)
          }
        }
      }
      j = j + 1
    }
    while (!is.null(dev.list()))
      dev.off()
  }
  )

}else{
  print("截取一个位点")
  localPos <- as.numeric(localPos)
  setMethod("chromatogram", "sangerseq", 
            function(obj, trim5=0, trim3=0, 
                    showcalls=c("primary", "secondary", "both", "none"), 
                    width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
                    filename=NULL, showtrim=FALSE, showhets=TRUE) {
              
    if (localPos + downStream_length > nchar(seq@primarySeq)){
      stop("截图长度不能超过原图总长度！")
    }
    if (localPos < upStream_length){
      stop("上游序列长度不能超过截取位置！")
    }
    
    upStream_start <- localPos - upStream_length
    upStream_end <- localPos -1
    downStream_start <- localPos + 1
    downStream_end <- localPos + downStream_length
    
    
    upStream_seq <- substr(seq@primarySeq,upStream_start,upStream_end)
    downStream_seq <- substr(seq@primarySeq,downStream_start,downStream_end)
    
    pattern <- paste0(upStream_seq,'([ACTG]*)',downStream_seq)
    
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    standSeq <- paste(basecalls1,collapse="")
    match <- str_extract_all(standSeq,pattern)
    
    
    revCom_up <- DNAStringSet(upStream_seq,)%>% reverse %>% complement
    revCom_down <- DNAStringSet(downStream_seq,)%>% reverse %>% complement
    revCom_pattern <- paste0(revCom_down,'([ACTG]*)',revCom_up)
    revCom_match  <- str_extract_all(standSeq,revCom_pattern)
    
    #target_seq <- str_extract_all(basecalls1,)
    #pic_str <- "CCTCCGTAAATACTGGACCCAAGTTACTGC" # wait_for_change
    pic_pos <- nchar(upStream_seq)+1 # wait_for_change
    #print(pic_pos)
    originalpar <- par(no.readonly=TRUE)
    #showcalls <- showcalls[3]
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    list_str_basecalls1 <- paste(basecalls1,collapse="")
    list_tmp_store_1 <- strsplit(list_str_basecalls1, match)
    list_tmp_store_2 <- strsplit(list_str_basecalls1, revCom_match)
    
    
    if (length(list_tmp_store_1[[1]]) == 2)
    {
      trim5 = nchar(list_tmp_store_1[[1]][1]) 
      trim3 = nchar(list_tmp_store_1[[1]][2])
    }
    else if (length(list_tmp_store_2[[1]]) == 2)
    {
      trim5 = nchar(list_tmp_store_2[[1]][1])
      trim3 = nchar(list_tmp_store_2[[1]][2])
      print(trim5)
      print(trim3)
      #pic_pos <- nchar(pic_str) - pic_pos + 1
    }
    aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(aveposition)] 
    basecalls2 <- basecalls2[1:length(aveposition)] 
    if(showtrim == FALSE) {
      if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
      else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
      if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
      else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
      aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
    }
    indexes <- 1:length(basecalls1)
    trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
    #false if not trimmed
    if (!is.null(trim3)) {
      traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                              nrow(traces))), ]
    }
    if (!is.null(trim5)) {
      offset <- max(c(1, aveposition[1] - 10))
      traces <- traces[offset:nrow(traces),]
      aveposition <- aveposition - (offset-1)
    }
    maxsignal <- apply(traces, 1, max)
    ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
    p <- c(0, aveposition, nrow(traces))
    midp <- diff(p)/2
    starts <- aveposition - midp[1:(length(midp)-1)]
    starthets <- starts
    #print(starthets)
    for(ii in 1:length(starthets))
    {
      if(ii != pic_pos)
      {
        starthets[ii] <- NA
      }
    }
    ends <- aveposition + midp[2:(length(midp))]
    endhets <- ends
    for(ii in 1:length(endhets))
    {
      if(ii != pic_pos)
      {
        endhets[ii] <- NA
      }
    }
    starttrims <- starts
    starttrims[!trimmed] <- NA
    endtrims <- ends
    endtrims[!trimmed] <- NA
    
    colortranslate <- c(A="green", C="blue", G="black", T="red")
    colorvector1 <- unname(colortranslate[basecalls1])
    colorvector1[is.na(colorvector1)] <- "purple"
    colorvector2 <- unname(colortranslate[basecalls2])
    colorvector2[is.na(colorvector2)] <- "purple"
    
    valuesperbase <- nrow(traces)/length(basecalls1)
    tracewidth <- width*valuesperbase
    breaks <- seq(1,nrow(traces), by=tracewidth)
    numplots <- length(breaks)
    if(!is.null(filename)) png(filename, width=2500, height=600, pointsize = 8, res=500) 
    par(mar=c(0,0,1,0), mfrow=c(numplots, 1))
    basecallwarning1 = 0
    basecallwarning2 = 0
    j = 1
    
    for(i in breaks) {
      range <- aveposition >= i & aveposition < (i+tracewidth)
      starthet <- starthets[range] - tracewidth*(j-1)
      starthet[starthet < 0] <- 0
      endhet <- endhets[range] - tracewidth*(j-1)
      endhet[endhet > tracewidth] <- tracewidth
      lab1 <- basecalls1[range]
      lab2 <- basecalls2[range]
      pos <- aveposition[range] - tracewidth*(j-1)
      colors1 <- colorvector1[range]
      colors2 <- colorvector2[range]
      starttrim <- starttrims[range] - tracewidth*(j-1)
      endtrim <- endtrims[range] - tracewidth*(j-1)
      plotrange <- i:min(i+tracewidth, nrow(traces))
      plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", bty="n", xlab="", yaxt="n", xlim=c(23,tracewidth-23))
      #title(main = "Stopping Distance versus Speed", font.main= 4) # wait_for_change
      if (showhets==TRUE) {
      #if (showhets==FALSE) {
        rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')#lty=3
      }
      if (showtrim==TRUE) {
        rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent', 
            density=15)
      }
      lines(traces[plotrange,1], col="green")
      lines(traces[plotrange,2], col="blue")
      lines(traces[plotrange,3], col="black")
      lines(traces[plotrange,4], col="red")
      #mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
      
      for(k in 1:length(lab1)) {
        if (showcalls=="primary" | showcalls=="both") {
          if (is.na(basecalls1[1]) & basecallwarning1==0) {
            warning("Primary basecalls missing")
            basecallwarning1 = 1
          } 
          else if (length(lab1) > 0) {   
            axis(side=3, at=pos[k], labels=lab1[k], cex.axis=1, col.axis=colors1[k], 
                family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0, 
                                                          -1), tick=FALSE)
          }
        }
        if (showcalls=="secondary" | showcalls=="both") {
          if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
            warning("Secondary basecalls missing")
            basecallwarning2 = 1
          } 
          else if (length(lab2) > 0) { 
            axis(side=3, at=pos[k], labels=lab2[k], cex.axis=1, col.axis=colors2[k], 
                family="mono", cex=cex.base, line=-1.0, tick=FALSE)
          }
        }
      }
      j = j + 1
    }
    while (!is.null(dev.list()))
      dev.off()
  }
  )
}
chromatogram(seq, trim5=0, trim3=0, showcalls="primary", width=41, height=4, cex.mtext=1, cex.base=1, ylim=3, filename=pngName, showtrim=FALSE, showhets=TRUE)
#chromatogram(seq, trim5=0, trim3=0, showcalls="both", width=70, height=4, cex.mtext=1, cex.base=1, ylim=3, filename=pngName, showtrim=FALSE, showhets=TRUE)
