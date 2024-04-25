#!/usr/bin/env Rscript
args=commandArgs(T) 
library(sangerseqR)
library(Biostrings)
library(stringr)
#library(tidyverse)


#args[1] == ab1 name; args[2] == png name; args[3] == upStreamSeq ; args[4] == downStreamSeq
#print(args[1])
seq <- readsangerseq(args[1])
#' @rdname chromatogram
setMethod("chromatogram", "sangerseq", 
  function(obj, trim5=0, trim3=0, 
           showcalls=c("primary", "secondary", "both", "none"), 
           width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
           filename=NULL, showtrim=FALSE, showhets=TRUE) {
    #
    upStream_seq <- args[3]
    downStream_seq <-  args[4]
    #upstream

    pattern_1 <- paste0(upStream_seq,'([ACTG]*)')
    #downstream
    pattern_2 <- paste0('([ACTG]*)',downStream_seq)
    
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    list_str_basecalls1 <- paste(basecalls1,collapse="")
    standSeq <- paste(basecalls1,collapse="")
    #(standSeq)
    #print(str_extract_all(standSeq, pattern_1, simplify = TRUE))
    #print(str_extract_all(standSeq, pattern_2, simplify = TRUE))
    if(length(str_extract_all(standSeq, pattern_1, simplify = TRUE)) != 0)
    {
      # match success
      match_tmp <- str_extract_all(standSeq,pattern_1)
      match <- str_sub(match_tmp, 1, 41)
      list_tmp_store_1 <- strsplit(list_str_basecalls1, match)
    }
    else if(length(str_extract_all(standSeq, pattern_2, simplify = TRUE)) != 0)
    {
      # match success
      match_tmp <- str_extract_all(standSeq,pattern_2)
      match <- str_sub(match_tmp, -41, -1)
      list_tmp_store_1 <- strsplit(list_str_basecalls1, match)
    }

    
    #
    revCom_up <- DNAStringSet(upStream_seq,)%>% reverse %>% complement
    revCom_down <- DNAStringSet(downStream_seq,)%>% reverse %>% complement
    #revCom_pattern <- paste0(revCom_down,'([ACTG]*)',revCom_up)
    #revCom_match  <- str_extract_all(standSeq,revCom_pattern)
    #upstream
    revCom_pattern_1 <- paste0('([ACTG]*)', revCom_up)
    #downstream
    revCom_pattern_2 <- paste0(revCom_down, '([ACTG]*)')
    
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    standSeq <- paste(basecalls1,collapse="")
    #print(str_extract_all(standSeq, revCom_pattern_1, simplify = TRUE))
    #print(str_extract_all(standSeq, revCom_pattern_2, simplify = TRUE))
    if(length(str_extract_all(standSeq, revCom_pattern_1, simplify = TRUE)) != 0)
    {
      # match success
      revMatch_tmp <- str_extract_all(standSeq,revCom_pattern_1)
      revCom_match <- str_sub(revMatch_tmp, -41, -1)
      list_tmp_store_2 <- strsplit(list_str_basecalls1, revCom_match)
    }
    else if(length(str_extract_all(standSeq, revCom_pattern_2, simplify = TRUE)) != 0)
    {
      # match success
      revMatch_tmp <- str_extract_all(standSeq,revCom_pattern_2)
      revCom_match <- str_sub(revMatch_tmp , 1, 41)
      list_tmp_store_2 <- strsplit(list_str_basecalls1, revCom_match)
    }

    #target_seq <- str_extract_all(basecalls1,)
    #pic_str <- "CCTCCGTAAATACTGGACCCAAGTTACTGC" # wait_for_change
    pic_pos <- 21 # wait_for_change
    originalpar <- par(no.readonly=TRUE)
    showcalls <- showcalls[1]
    traces <- obj@traceMatrix
    #basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    #basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    #list_str_basecalls1 <- paste(basecalls1,collapse="")
    #list_tmp_store_1 <- strsplit(list_str_basecalls1, match)
    #list_tmp_store_2 <- strsplit(list_str_basecalls1, revCom_match)
    #if (length(list_tmp_store_1[[1]]) == 2)
    if (exists("list_tmp_store_1"))
    {
      trim5 = nchar(list_tmp_store_1[[1]][1]) 
      trim3 = nchar(list_tmp_store_1[[1]][2])
    }
    #else if (length(list_tmp_store_2[[1]]) == 2)
    else if((exists("list_tmp_store_2")))
    {
      trim5 = nchar(list_tmp_store_2[[1]][1])
      trim3 = nchar(list_tmp_store_2[[1]][2])
      #print(trim5)
      #print(trim3)
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
    #par(mar=c(0,0,1,0), mfrow=c(numplots, 1))
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
chromatogram(seq, trim5=0, trim3=0, showcalls=c("primary", "secondary", "both", "none"), width=41, height=4, cex.mtext=1, cex.base=1, ylim=3, filename=args[2], showtrim=FALSE, showhets=TRUE)
