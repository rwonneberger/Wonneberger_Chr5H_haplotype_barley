###########################################################################
## Dotplots aligning the chromosome 5H region between Barke and Morex V2 ##
###########################################################################

#This script plots pairwise dotplots. It uses a modified wrapper script for lastz written by Mark Timothy Rabanus-Wallace, IPK (now University of Melbourne)

source("https://raw.githubusercontent.com/mtrw/tim_r_functions/master/tim_functions.R")
source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

get_lastz_dotplot2_allaxes <- function(
  file1,
  file2,
  name1,
  name2,
  range1=NULL,
  range2=NULL,
  seq1=NULL,
  seq2=NULL,
  annot1=NULL,
  annot2=NULL,
  lastz_binary=""/path/to/your/lastz/"",
  min_length_plot=0,
  save_alignments_to_file=NULL,
  save_dots_to_file=NULL,
  plot_from_file=NULL,
  args="--notransition --step=150 --nogapped",
  plot = F
){
  
  if(is.null(plot_from_file)){ #we need to make the alignments
    
    tfo <- tempfile() #output (alignment)
    tfd <- tempfile() #output (dotplot info)
    
    file1call <- file1
    file2call <- file2
    
    options(scipen = 999)
    if(!is.null(seq1)){
      tf1 <- tempfile()
      writeLines(seq1,tf1)
      file1call <- paste0(file1call,"[subset=",tf1)
      if(!is.null(range1)){
        file1call <- paste0(file1call,",",range1[1],"..",range1[2])
      }
      file1call <- paste0(file1call,"]")
    }
    if(!is.null(seq2)){
      tf2 <- tempfile()
      writeLines(seq2,tf2)
      file2call <- paste0(file2call,"[subset=",tf2)
      if(!is.null(range2)){
        file2call <- paste0(file2call,",",range2[1],"..",range2[2])
      }
      file2call <- paste0(file2call,"]")
    }
    
    options(scipen = 0)
    cmd <- paste0(lastz_binary," ",file1call," ",file2call," ",args," --rdotplot=",tfd," > ",tfo)
    #system(paste("cat",tf1))
    #system(paste("cat",tf2))
    ce("Running command: ",cmd)
    system(cmd)
    
    if(!is.null(save_alignments_to_file)){
      file.copy(tfo,save_alignments_to_file)
      ce("Alignments saved as ",save_alignments_to_file," in ",getwd())
    }
    
    if(!is.null(save_dots_to_file)){
      file.copy(tfd,save_dots_to_file)
      ce("Dots saved as ",save_dots_to_file," in ",getwd())
    }
    
    unlink(tf1)
    unlink(tf2)
    unlink(tfo)
    dp <- fread(tfd,header=T,col.names=c("s1","s2"))
    unlink(tfd)
  } else {
    tfd <- plot_from_file
    dp <- fread(tfd,header=T,col.names=c("s1","s2"))
  }
  
  suppressWarnings(dp[,s1:=as.numeric(s1)])
  suppressWarnings(dp[,s2:=as.numeric(s2)])
  dp[,idx:=(1:.N)%%3]
  #dev dp <- dp[1:54]
  abs(dp[idx==2,s1]-dp[idx==1,s1]) -> l1_
  abs(dp[idx==2,s2]-dp[idx==1,s2]) -> l2_
  dp[,l1:=rep(l1_,each=3)]
  dp[,l2:=rep(l2_,each=3)]
  dp[,l:=pmax(l1,l2)]
  dp <- dp[l>=min_length_plot]
  dp[l!=max(l) & idx!=0,c:=replace_scale_with_colours(-log(l))] #,fun="sequential_hcl",palette="Reds 3"
  dp[l==max(l) & idx!=0,c:="#000000"]
  
  seq1descript <- paste0(name1) 
  seq2descript <- paste0(name2) 

  
  if(plot==F){
    return()
  }
  #dev.off()
  #par(mar=c(0,0,0,0))
  par(mar=c(4,4,1,1))
  
  null_plot(
    x=dp$s1,
    y=dp$s2,
    xaxt='n',
    yaxt='n', 
    xlab=seq1descript,
    ylab=seq2descript 
  )
  x0=dp$s1
  y0=dp$s2
  axis(1,at=axTicks(1), labels=sprintf("%.0f", axTicks(1)/1000000)) 
  axis(2,at=axTicks(2), labels=sprintf("%.0f", axTicks(2)/1000000))
  
  l_ply(seq(from=1,length.out=nrow(dp)/3,by=3),function(i){
    lines(
      x=dp[i:(i+2),s1],
      y=dp[i:(i+2),s2],
      col=dp[i:(i+2),c],
      lwd=1
    )
  })
  
  if(!is.null(annot1)){
    fannot1 <- annot1[seqname==seq1 & feature=="exon" & ((end %between% range(dp$s1,na.rm=T)) | (start %between% range(dp$s1,na.rm=T))) ]
    
    if (nrow(fannot1)>0){
      
      if(is.null(annot1$col)){
        fannot1[,col:="#f02222"]
      }
      if(is.null(annot1$linecol)){
        fannot1[,linecol:="#f0222211"]
      }
      
      fannot1[,idx:=1:.N]
      pannot1 <- fannot1[,{
        data.table(
          x=c(start,end,NA),
          y=min(dp$s2,na.rm=T)+(0.02)*abs(diff(range(dp$s2,na.rm=T))),
          c=col,
          lc=linecol
        )
      },by="idx"]
      pannot1[is.na(x),y:=NA][is.na(y),c:=NA][is.na(c),lc:=NA]
      
      l_ply(seq(from=1,length.out=nrow(pannot1)/3,by=3),function(i){
        lines(
          x=pannot1[i:(i+2),x],
          y=pannot1[i:(i+2),y],
          col=pannot1[i:(i+2),c],
          lwd=4,
          lend="butt"
        )
      })
      
      abline(
        v=pannot1[!is.na(x),x],
        col=pannot1[!is.na(x),lc],
        lwd=.5
      )
    }
  }
  
  if(!is.null(annot2)){
    fannot2 <- annot2[seqname==seq2 & feature=="exon" & ((end %between% range(dp$s2,na.rm=T)) | (start %between% range(dp$s2,na.rm=T))) ]
    
    if (nrow(fannot2)>0){
      
      if(is.null(annot2$col)){
        fannot2[,col:="#f02222"]
      }
      if(is.null(annot2$linecol)){
        fannot2[,linecol:="#f0222211"]
      }
      
      fannot2[,idx:=1:.N]
      pannot2 <- fannot2[,{
        data.table(
          y=c(start,end,NA),
          x=min(dp$s1,na.rm=T)+(0.02)*abs(diff(range(dp$s1,na.rm=T))),
          c=col,
          lc=linecol
        )
      },by="idx"]
      pannot2[is.na(x),y:=NA][is.na(y),c:=NA][is.na(c),lc:=NA]
      
      l_ply(seq(from=1,length.out=nrow(pannot2)/3,by=3),function(i){
        lines( #bars
          x=pannot2[i:(i+2),x],
          y=pannot2[i:(i+2),y],
          col=pannot2[i:(i+2),c],
          lwd=4,
          lend="butt"
        )
      })
      
      abline( #lines
        h=pannot2[!is.na(x),y],
        col=pannot2[!is.na(x),lc],
        lwd=.5
      )
    }
  }
  
  
}


# Note: this script makes the alignments first and saves them to file, then makes the figures. You can do this in one step
################ paths #########################

#####################
## Make alignments ##
#####################

generalpath='/add/your/path/here/'

pathMorex=paste0(generalpath, 'Morex_pseudomolecules_v2.fasta')
pathBarke=paste0(generalpath, 'Barke_pseudomolecules_v1.fasta')


rangeMorex=c(46242159,325574358)
rangeBarke=c(48000000,330000000)


#Global parameters
seq1="chr5H"
seq2="chr5H"
min_length_plot=10000 #5000
annot1=NULL
annot2=NULL
args = "--notransition --step=100000 --nogapped" #10000
lastz_binary="/path/to/your/lastz/"


file1= pathBarke
name1="Barke"
range1 = rangeBarke

file2=pathMorex
name2="Morex"
range2=rangeMorex
get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2, name1 = name1,  name2 = name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_5H_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_5H_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary, plot = F)



#Global params
seq1="chr5H"
seq2="chr5H"
min_length_plot=1000 #1000
annot1=NULL
annot2=NULL
args = "--notransition --step=4000 --nogapped" #4000
lastz_binary="/path/to/your/lastz/"

rangeMorex=c(55000000,82000000)
rangeBarke=c(54000000,86000000)


file1= pathBarke
name1="Barke"
range1 = rangeBarke

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2, name1 = name1,  name2 = name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_5H_start_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_5H_start_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary, plot = F)


#Global params
seq1="chr5H"
seq2="chr5H"
min_length_plot=500 #500
annot1=NULL
annot2=NULL
args = "--notransition --step=1000 --nogapped" # 1000
lastz_binary="/path/to/your/lastz/"

rangeMorex=c(305000000,318000000)
rangeBarke=c(307000000,322000000)


file1= pathBarke
name1="Barke"
range1 = rangeBarke

file2=pathMorex
name2="Morex"
range2=rangeMorex


get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2, name1 = name1,  name2 = name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_5H_end_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_5H_end_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary, plot = F)



#################
## Make figure ##
#################

generalpath='/add/your/path/here/'

pathMorex=paste0(generalpath, 'Morex_pseudomolecules_v2.fasta')
pathBarke=paste0(generalpath, 'Barke_pseudomolecules_v1.fasta')

rangeMorex=c(46242159,325574358)
rangeBarke=c(48000000,330000000)


#Global params
seq1="chr5H"
seq2="chr5H"
min_length_plot=10000 #5000
annot1=NULL
annot2=NULL
args = "--notransition --step=300000 --nogapped" #10000
lastz_binary="/path/to/your/lastz/"


file1= pathBarke
name1="Barke"
range1 = rangeBarke

file2=pathMorex
name2="Morex"
range2=rangeMorex

plot1<-ggdraw(function () get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_5H_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary, plot = T))


#Global params
seq1="chr5H"
seq2="chr5H"
min_length_plot=2000 #1000
annot1=NULL
annot2=NULL
args = "--notransition --step=5000 --nogapped" #4000
lastz_binary="/path/to/your/lastz/"

rangeMorex=c(55000000,82000000)
rangeBarke=c(54000000,86000000)


file1= pathBarke
name1="Barke"
range1 = rangeBarke

file2=pathMorex
name2="Morex"
range2=rangeMorex


plot2<-ggdraw(function () get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_5H_start_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary, plot = T))


#Global params
seq1="chr5H"
seq2="chr5H"
min_length_plot=1000 #500
annot1=NULL
annot2=NULL
args = "--notransition --step=2000 --nogapped" # 1000
lastz_binary="/path/to/your/lastz/"


rangeMorex=c(305000000,318000000)
rangeBarke=c(307000000,322000000)


file1= pathBarke
name1="Barke"
range1 = rangeBarke

file2=pathMorex
name2="Morex"
range2=rangeMorex



plot3<-ggdraw(function () get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_5H_end_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary, plot = T))


tiff("ESM18.tiff", width=45, height=15, res=300, units="cm")

ggarrange(plot1, plot2, plot3, ncol = 3, labels = c("a", "b", "c"))
dev.off()
