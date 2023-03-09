###################################################################################################
## Dotplots aligning the chromosome 5H region between Barke and the accessions from pangenome V1 ##
###################################################################################################

#This script plots pairwise dotplots. It uses a modified wrapper script for lastz written by Mark Timothy Rabanus-Wallace, IPK (now University of Melbourne)

source("https://raw.githubusercontent.com/mtrw/tim_r_functions/master/tim_functions.R")
source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

################ paths #########################
rangeMorex=c(50949849,325574358)
rangeGoldenPromise=c(41696520,298782928)
rangeHOR13821=c(55478598,333876948)
rangeHOR3081=c(51219855,324538253)
rangeHOR3365=c(54696442,344810825)

rangeBarke=c(50000000,330000000)
rangeIgri=c(50461764,328340288)
rangeHOR13942=c(51590962,330933960)
rangePlanet=c(50177890,336962473)

generalpath='/add/your/path/here/'

pathMorex=paste0(generalpath, 'Morex_pseudomolecules_v2.fasta')
pathGoldenPromise=paste0(generalpath, 'GoldenPromise_pseudomolecule_v1.fasta')
pathHOR13821=paste0(generalpath, 'HOR13821_pseudomolecule_v1.fasta')
pathHOR3081=paste0(generalpath, 'HOR3081_pseudomolecule_v1.fasta')
pathHOR3365=paste0(generalpath, 'HOR3365_pseudomolecule_v1.fasta')

pathBarke=paste0(generalpath, 'Barke_pseudomolecules_v1.fasta')
pathIgri=paste0(generalpath, 'Igri_pseudomolecule_v1.fasta')
pathHOR13942=paste0(generalpath, 'HOR13942_pseudomolecule_v1.fasta')
pathPlanet =paste0(generalpath, 'RGTPlanet_pseudomolecule_v1.fasta')


get_lastz_dotplot <- function(
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
  lastz_binary="/path/to/your/lastz/",
  min_length_plot=0,
  save_alignments_to_file=NULL,
  save_dots_to_file=NULL,
  plot_from_file=NULL,
  args="--notransition --step=150 --nogapped",
  plotname="dummy"
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
  tiff(paste0(plotname, ".tiff"), height=20, width=20, units="cm", res=300)
  
  #dev.off()
  #  par(mar=c(0,0,0,0))
  par(mar=c(5,5,2,2))
  
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
  dev.off()
  
}


################ paths #########################
rangeMorex=c(50949849,325574358)
rangeGoldenPromise=c(41696520,298782928)
rangeHOR13821=c(55478598,333876948)
rangeHOR3081=c(51219855,324538253)
rangeHOR3365=c(54696442,344810825)

rangeBarke=c(50000000,330000000)
rangeIgri=c(50461764,328340288)
rangeHOR13942=c(51590962,330933960)
rangePlanet=c(50177890,336962473)


generalpath='/add/your/path/here/'

pathMorex=paste0(generalpath, 'Morex_pseudomolecules_v2.fasta')
pathGoldenPromise=paste0(generalpath, 'GoldenPromise_pseudomolecule_v1.fasta')
pathHOR13821=paste0(generalpath, 'HOR13821_pseudomolecule_v1.fasta')
pathHOR3081=paste0(generalpath, 'HOR3081_pseudomolecule_v1.fasta')
pathHOR3365=paste0(generalpath, 'HOR3365_pseudomolecule_v1.fasta')

pathBarke=paste0(generalpath, 'Barke_pseudomolecules_v1.fasta')
pathIgri=paste0(generalpath, 'Igri_pseudomolecule_v1.fasta')
pathHOR13942=paste0(generalpath, 'HOR13942_pseudomolecule_v1.fasta')
pathPlanet =paste0(generalpath, 'RGTPlanet_pseudomolecule_v1.fasta')


#Global params
seq1="chr5H"
seq2="chr5H" 
min_length_plot=1000 #1000
annot1=NULL
annot2=NULL
args = "--notransition --step=50000 --nogapped" #50000
lastz_binary="/path/to/your/lastz/"




# GP 
file1= pathGoldenPromise
name1="GoldenPromise"

range1 = rangeGoldenPromise

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))




# Morex
file1= pathMorex
name1="Morex"

range1 = rangeMorex

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))



# HOR13821
file1= pathHOR13821
name1="HOR13821"

range1 = rangeHOR13821

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))





# HOR3081
file1= pathHOR3081
name1="HOR3081"

range1 = rangeHOR3081

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))





# HOR3365
file1= pathHOR3365
name1="HOR3365"

range1 = rangeHOR3365

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))





# Barke
file1= pathBarke
name1="Barke"

range1 = rangeBarke

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))





# Igri
file1= pathIgri
name1="Igri"

range1 = rangeIgri

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))





# HOR13942
file1= pathHOR13942
name1="HOR13942"

range1 = rangeHOR13942

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))





# Planet
file1= pathPlanet
name1="RGTPlanet"

range1 = rangePlanet

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, plot_from_file=NULL,  save_alignments_to_file = paste0("Results/", name1, "_", name2, "_align"),  save_dots_to_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary,  plotname=paste0("Figs/", name1, "_", name2))







#Global params
seq1="chr5H"
seq2="chr5H" 
min_length_plot=2000 #1000
annot1=NULL
annot2=NULL
args = "--notransition --step=60000 --nogapped" #50000
lastz_binary="/path/to/lastz"






## This fct plots both axes

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
  lastz_binary="/path/to/your/lastz/",
  min_length_plot=0,
  save_alignments_to_file=NULL,
  save_dots_to_file=NULL,
  plot_from_file=NULL,
  args="--notransition --step=150 --nogapped"
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
  
  par(mar=c(4,4,1,1))
  
  null_plot(
    x=dp$s1,
    y=dp$s2,
    xaxt='n',

    yaxt='n', 
    xlab=seq1descript,
    ylab=seq2descript,
    cex.lab=1.3 
  )
  x0=dp$s1
  y0=dp$s2
  axis(1,at=axTicks(1), labels=sprintf("%.1f", axTicks(1)/1000000)) 
  axis(2,at=axTicks(2), labels=sprintf("%.1f", axTicks(2)/1000000))
  
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



## This fct plots only the x axis

get_lastz_dotplot2_x <- function(
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
  lastz_binary="/path/to/your/lastz/",
  min_length_plot=0,
  save_alignments_to_file=NULL,
  save_dots_to_file=NULL,
  plot_from_file=NULL,
  args="--notransition --step=150 --nogapped"
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
  
  par(mar=c(4,4,1,1))
  
  null_plot(
    x=dp$s1,
    y=dp$s2,
    xaxt='n', 

    yaxt='n',
    xlab=seq1descript,
    ylab=' ',
    cex.lab=1.3 
  )
  x0=dp$s1
  #y0=dp$s2
  axis(1,at=axTicks(1), labels=sprintf("%.1f", axTicks(1)/1000000))

   
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






## This fct plots only the y axis

get_lastz_dotplot2_y <- function(
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
  lastz_binary="/path/to/your/lastz/",
  min_length_plot=0,
  save_alignments_to_file=NULL,
  save_dots_to_file=NULL,
  plot_from_file=NULL,
  args="--notransition --step=150 --nogapped"
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
  
  seq2descript <- paste0(name2) 
  
  
  par(mar=c(4,4,1,1))
  
  null_plot(
    x=dp$s1,
    y=dp$s2,
    xaxt='n', 
    yaxt='n', 
    xlab=' ',
    ylab=seq2descript,
    cex.lab=1.3 
  )
  #x0=dp$s1
  y0=dp$s2
  
  axis(2,at=axTicks(2), labels=sprintf("%.1f", axTicks(2)/1000000)) 
  
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






## This fct doesn't plot any axes

get_lastz_dotplot2_none <- function(
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
  lastz_binary="/path/to/your/lastz/",
  min_length_plot=0,
  save_alignments_to_file=NULL,
  save_dots_to_file=NULL,
  plot_from_file=NULL,
  args="--notransition --step=150 --nogapped"
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
  

  
  
  par(mar=c(4,4,1,1))
  
  null_plot(
    x=dp$s1,
    y=dp$s2,
    xaxt='n',
    yaxt='n', 
    xlab=' ',
    ylab=' ',
    cex.lab=1.3 
    
  )
  x0=dp$s1
  y0=dp$s2

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



print("Start plotting")

tiff("Figs/ESM19.tiff", width=80, height=80, res=300, units="cm")

par(mfcol=c(9,9))


# GP 
file1= pathGoldenPromise
name1="GoldenPromise"

range1 = rangeGoldenPromise

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_y(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_allaxes(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))



print("GP done")




# Morex
file1= pathMorex
name1="Morex"

range1 = rangeMorex


file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

print("Morex done")

# HOR13821
file1= pathHOR13821
name1="HOR13821"

range1 = rangeHOR13821

range1 = rangeMorex

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

print("13821 done")

# HOR3081
file1= pathHOR3081
name1="HOR3081"

range1 = rangeHOR3081

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))


print("3081 done")

# HOR3365
file1= pathHOR3365
name1="HOR3365"

range1 = rangeHOR3365

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

print("3365 done")

# Barke
file1= pathBarke
name1="Barke"

range1 = rangeBarke

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

print("Barke done")

# Igri
file1= pathIgri
name1="Igri"

range1 = rangeIgri

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

print("Igii done")

# HOR13942
file1= pathHOR13942
name1="HOR13942"

range1 = rangeHOR13942

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))



print("13942 done")



# Planet
file1= pathPlanet
name1="RGTPlanet"

range1 = rangePlanet

file2=pathGoldenPromise
name2="GoldenPromise"
range2=rangeGoldenPromise

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathMorex
name2="Morex"
range2=rangeMorex

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13821
name2="HOR13821"
range2=rangeHOR13821

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3081
name2="HOR3081"
range2=rangeHOR3081

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR3365
name2="HOR3365"
range2=rangeHOR3365

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathBarke
name2="Barke"
range2=rangeBarke

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathIgri
name2="Igri"
range2=rangeIgri

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathHOR13942
name2="HOR13942"
range2=rangeHOR13942

get_lastz_dotplot2_none(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

file2=pathPlanet
name2="RGTPlanet"
range2=rangePlanet

get_lastz_dotplot2_x(  file1 = file1,  file2 = file2,  name1= name1,  name2= name2,  seq1 = seq1,  seq2 = seq2,  range1 = range1,  range2 = range2,  args = args, save_dots_to_file=NULL,  save_alignments_to_file = NULL,  plot_from_file = paste0("Results/", name1, "_", name2, "_dot"), min_length_plot=min_length_plot,  annot1 = annot1, annot2 = annot2,  lastz_binary=lastz_binary)

print(paste0("Done with ", name1, " vs ", name2))

print("Planet done")
dev.off()

print("all done")

