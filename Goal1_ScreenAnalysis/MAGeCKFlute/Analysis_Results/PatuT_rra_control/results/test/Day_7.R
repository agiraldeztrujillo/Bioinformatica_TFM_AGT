pdf(file='Day_7.pdf',width=4.5,height=4.5);
gstable=read.table('Day_7.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("PADI2","SMARCA2","PRMT3","TDRKH","BPTF","IWS1","HELLS","ING4","CHD3","ASXL1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='day7_Cas9_r1,day7_Cas9_r2,day7_Cas9_r3,day7_Cas9_r4,day7_Cas9_r5_vs_day7_r1,day7_r2,day7_r3,day7_r4,day7_r5 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(103.01098379659174,77.92692042368311,99.7760680888569,73.65473721572297,99.21810564977768,105.1111633218437,90.6048065586785,105.82961777883496,74.76002406543209,89.42417136233267),c(61.44514822954595,49.808959446065494,73.93478426728245,49.694762458801044,60.07701809986539,33.486742297224545,47.39328343069337,57.7252460611827,77.02547934014216,59.08382750725551),c(96.6857479494326,56.23592195523524,77.52385146472335,87.85324077538041,58.25650239986947,40.92824058549667,110.11968797131695,67.34612040471316,27.18546329652076,31.937204057975954),c(79.51725065000065,89.97747512837638,94.75137401243964,72.76733074324439,96.48733209978381,85.5772303151294,131.02848948485814,109.43744565765888,101.94548736195286,105.39277339132065))
targetgene="PADI2"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(83.13167113409159,58.64603289617389,47.37568700621982,86.96583430290183,86.47449574980624,49.2999261598028,66.90816484333182,55.320027475300094,55.50365423039655,75.0524295362435),c(63.25235847159142,49.808959446065494,61.7319557959834,43.48291715145091,78.2821750998246,59.53198630617697,66.90816484333182,75.7643854553023,77.02547934014216,87.82731115943388),c(80.42085577102338,84.35388293285285,58.86070203803069,58.56882718358695,52.79495529988171,62.32254816427901,33.45408242166591,58.92785535412401,70.22911351601196,76.64928973914229),c(56.02351750340954,36.15166411407979,47.37568700621982,66.55548543589425,64.62830734985519,36.27730415532659,30.666242219860415,58.92785535412401,77.02547934014216,27.14662344927956))
targetgene="SMARCA2"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(64.15596359261416,69.08984697357472,63.88539611444794,93.17767961025196,75.55140154983071,90.22816674529946,71.08992514604006,68.54872969765447,52.10547131833146,55.89010710145792),c(65.05956871363689,61.05614383711254,61.01414235649522,57.68142071110835,61.89753379986131,37.2074914413606,26.484481917152177,46.901762424710945,46.4418331315563,47.905806086963935),c(5.421630726136407,4.820221881877306,2.8712537579527164,7.0992517798287205,0.9102578499979604,0.0,8.363520605416477,1.2026092929413064,3.398182912065095,0.0),c(36.14420484090938,34.54492348678736,36.608485413897135,34.60885242666501,30.948766899930654,35.34711686929258,16.727041210832954,52.91480888941748,23.787280384455666,22.356042840583168))
targetgene="PRMT3"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(46.9874662931822,85.15725324649908,67.47446331188884,67.44289190837284,73.73088584983479,47.43955158773477,30.666242219860415,48.104371717652256,67.9636582413019,73.4555693333447),c(46.083861172159466,54.6291813279428,50.246940764172535,44.3703236239295,59.16676024986743,69.76404645255113,50.18112363249887,45.69915313176964,46.4418331315563,41.51836527536874),c(59.637937987500486,48.20221881877306,61.01414235649522,61.23104660102271,73.73088584983479,30.6961804391225,43.21152312798513,38.483497374121804,27.18546329652076,31.937204057975954),c(0.0,0.0,0.0,0.0,4.551289249989802,1.8603745720680303,0.0,0.0,0.0,0.0))
targetgene="TDRKH"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(46.9874662931822,51.41570007335793,25.123470382086268,44.3703236239295,47.333408199893945,68.83385916651712,16.727041210832954,36.07827878823919,24.9200080218107,25.549763246380763),c(0.0,0.8033703136462177,0.0,1.7748129449571801,0.0,0.0,0.0,1.2026092929413064,0.0,0.0))
targetgene="BPTF"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(63.25235847159142,25.707850036678966,31.58379133747988,73.65473721572297,66.44882304985111,66.97348459444909,33.45408242166591,31.267841616473966,31.716373845940886,47.905806086963935),c(20.782917783522898,28.11796097761762,27.276910700550808,16.86072297709321,21.84618839995105,20.464120292748333,18.1209613117357,15.633920808236983,9.061821098840253,57.48696730435672),c(39.758625325000324,32.93818285949492,43.78661980877892,44.3703236239295,20.02567269995513,26.045244008952423,50.18112363249887,45.69915313176964,36.247284395361014,17.565462231886777),c(93.9749325863644,48.20221881877306,58.86070203803069,78.09176957811593,62.80779164985927,70.69423373858515,96.18048696228949,58.92785535412401,77.02547934014216,59.08382750725551))
targetgene="IWS1"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(120.17948109602371,98.01117826483856,119.87484439452591,149.08428737640313,119.24377834973282,93.01872860340151,101.75616736590047,84.18265050589144,79.29093461485222,172.46090191307016),c(47.89107141420494,41.77525630960332,47.37568700621982,33.72144595418642,36.410313999918415,26.045244008952423,23.696641715346686,36.07827878823919,23.787280384455666,62.27754791305311),c(70.4811994397733,40.9718859959571,38.761925732361675,39.045884789057965,30.038509049932696,47.43955158773477,69.69600504513731,43.293934545887026,43.043650219491205,59.08382750725551),c(35.24059971988665,30.52807191855627,38.04411229287349,39.045884789057965,28.217993349936773,39.99805329946265,26.484481917152177,30.06523232353266,45.30910549420127,39.92150507246994))
targetgene="HELLS"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(52.40909701931861,38.56177505501845,35.17285853492078,55.01920129367258,65.53856519985315,45.579177015666744,50.18112363249887,18.039139394119594,21.521825109745603,36.72778466667235),c(64.15596359261416,88.37073450108394,57.42507515905433,83.41620841298746,77.37191724982664,50.230113445836814,66.90816484333182,67.34612040471316,64.5654753292368,33.53406426087475),c(42.469440688068524,40.9718859959571,43.06880636929075,43.48291715145091,45.51289249989802,25.115056722918407,44.60544322888788,31.267841616473966,43.043650219491205,44.712085681166336),c(6.325235847159142,1.6067406272924354,2.153440318464537,5.324438834871541,0.0,0.9301872860340151,0.0,1.2026092929413064,0.0,0.0))
targetgene="ING4"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(60.54154310852322,88.37073450108394,53.11819452212525,109.15099611486657,82.8334643498144,59.53198630617697,36.2419226234714,90.19569697059798,80.42366225220725,55.89010710145792),c(47.89107141420494,41.77525630960332,57.42507515905433,55.906607766151176,55.525728849875584,37.2074914413606,51.57504373340161,48.104371717652256,38.51273967007108,25.549763246380763),c(55.11991238238681,42.57862662324954,50.246940764172535,55.906607766151176,36.410313999918415,43.71880244359871,37.63584272437415,38.483497374121804,60.03456477981668,41.51836527536874),c(48.794676535227666,35.34829380043358,52.400381082637075,25.73478770187911,23.666704099946973,40.92824058549667,32.06016232076316,30.06523232353266,39.64546730742611,35.13092446377355))
targetgene="CHD3"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(127.40832206420558,94.79769701025369,96.90481433090417,52.356981876236816,78.2821750998246,72.55460831065318,62.72640454062358,101.01918060706973,106.47639791137298,67.0681285217495),c(31.62617923579571,26.511220350325186,43.78661980877892,36.38366537162219,28.217993349936773,15.813183862578256,34.848002522568656,24.052185858826128,21.521825109745603,43.11522547826754),c(37.04780996193212,27.3145906639714,46.657873566731645,40.82069773401514,42.78211894990414,20.464120292748333,25.090561816249433,39.68610666706311,16.990914560325475,33.53406426087475),c(34.33699459886392,29.724701604910056,11.485015031810866,22.18516181196475,23.666704099946973,26.045244008952423,20.908801513541192,22.84957656588482,9.061821098840253,14.37174182608918))
targetgene="ASXL1"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("SETD6","SP140L","PRKAA1","PHF21B","NAP1L3","SETD7","PRDM10","DOT1L","NSD1","SP140")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='day7_Cas9_r1,day7_Cas9_r2,day7_Cas9_r3,day7_Cas9_r4,day7_Cas9_r5_vs_day7_r1,day7_r2,day7_r3,day7_r4,day7_r5 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(45.18025605113673,34.54492348678736,36.608485413897135,62.11845307350131,55.525728849875584,41.85842787153068,58.54464423791534,37.2808880811805,56.63638186775158,65.4712683188507),c(10.843261452272815,24.904479723032747,11.485015031810866,12.423690614700261,17.29489914996125,20.464120292748333,39.02976282527689,13.22870222235437,58.901837142461645,35.13092446377355))
targetgene="SETD6"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(60.54154310852322,37.75840474137223,43.78661980877892,33.72144595418642,57.34624454987151,60.46217359221098,34.848002522568656,27.660013737650047,32.84910148329592,54.29324689855912),c(64.15596359261416,47.39884850512684,68.91009019086519,47.03254304136527,50.06418174988782,79.06591931289128,78.05952565055378,63.738292525889236,35.11455675800598,71.8587091304459),c(49.6982816562504,30.52807191855627,54.553821401101615,70.10511132580861,50.974439599885784,63.252735450313025,30.666242219860415,39.68610666706311,61.16729241717171,59.08382750725551),c(33.43338947784118,30.52807191855627,41.63317949031439,27.50960064683629,28.217993349936773,57.671611734108936,44.60544322888788,31.267841616473966,54.37092659304152,59.08382750725551))
targetgene="SP140L"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(23.4937331465911,59.44940320982011,40.91536605082621,39.93329126153655,43.6923767999021,39.067866013428635,52.968963834304354,56.5226367682414,27.18546329652076,73.4555693333447),c(20.782917783522898,25.707850036678966,39.47973917184985,13.311097087178851,46.423150349895984,37.2074914413606,65.51424474242907,16.83653010117829,40.77819494478114,38.324644869571145),c(11.74686657329555,12.050554704693266,5.742507515905433,8.874064724785901,11.833352049973486,12.092434718442197,13.939201009027462,0.0,21.521825109745603,6.387440811595191),c(23.4937331465911,8.837073450108395,24.40565694259809,9.761471197264491,13.653867749969407,15.813183862578256,37.63584272437415,8.418265050589145,24.9200080218107,33.53406426087475))
targetgene="PRKAA1"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(70.4811994397733,78.73029073732934,94.75137401243964,55.01920129367258,54.61547099987763,60.46217359221098,52.968963834304354,54.11741818235878,46.4418331315563,63.87440811595191),c(58.73433286647775,66.67973603263607,49.52912732468436,40.82069773401514,54.61547099987763,51.160300731870834,52.968963834304354,27.660013737650047,37.380012032716046,62.27754791305311),c(26.204548509659304,27.3145906639714,22.97003006362173,25.73478770187911,20.935930549953092,38.13767872739462,19.514881412638445,13.22870222235437,20.38909747239057,12.774881623190382),c(16.264892178409223,3.2134812545848708,0.7178134394881791,2.6622194174357703,10.923094199975525,3.7207491441360605,16.727041210832954,22.84957656588482,7.9290934614852215,14.37174182608918))
targetgene="PHF21B"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(57.830727745455015,45.79210787783441,45.22224668775528,33.72144595418642,41.871861099906184,51.160300731870834,64.12032464152632,58.92785535412401,55.50365423039655,59.08382750725551),c(82.22806601306885,88.37073450108394,89.00886649653421,52.356981876236816,90.11552714979808,84.64704302909537,105.93792766860871,68.54872969765447,78.15820697749719,84.63359075363628),c(28.01175875170477,27.3145906639714,22.97003006362173,23.072568284443342,28.217993349936773,47.43955158773477,39.02976282527689,30.06523232353266,43.043650219491205,36.72778466667235),c(142.76960912159205,113.2752142241167,151.4586357320058,103.82655727999504,101.0386213497736,133.94696918889818,153.3312110993021,123.86875717295456,77.02547934014216,164.47660089857615))
targetgene="NAP1L3"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(95.78214282840987,61.85951415075876,90.44449337551056,130.44875145435273,81.01294864981848,124.64509632855803,71.08992514604006,64.94090181883054,60.03456477981668,67.0681285217495),c(115.66145549091003,126.93250955610239,91.88012025448693,86.96583430290183,91.02578499979604,106.97153789391174,85.02912615506752,104.62700848589365,104.21094263666292,79.84301014493988),c(61.44514822954595,69.08984697357472,33.01941821645624,25.73478770187911,39.1410875499123,80.92629388495932,75.2716854487483,91.39830626353928,40.77819494478114,41.51836527536874),c(5.421630726136407,10.44381407740083,0.0,0.0,0.0,0.0,1.3939201009027462,0.0,0.0,0.0))
targetgene="SETD7"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(8.132446089204612,20.88762815480166,17.945335987204476,11.536284142221671,15.474383449965327,29.765993153088484,16.727041210832954,27.660013737650047,27.18546329652076,4.790580608696393),c(18.07210242045469,30.52807191855627,39.47973917184985,21.297755339486162,26.397477649940853,18.6037457206803,9.757440706319223,13.22870222235437,28.31819093387579,28.74348365217836),c(57.830727745455015,40.9718859959571,78.9594783436997,33.72144595418642,43.6923767999021,59.53198630617697,41.817603027082384,25.254795151767432,54.37092659304152,30.34034385507716),c(18.975707541477426,8.837073450108395,14.356268789763583,23.95997475692193,15.474383449965327,26.97543129498644,18.1209613117357,24.052185858826128,27.18546329652076,27.14662344927956))
targetgene="PRDM10"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(45.18025605113673,32.93818285949492,33.01941821645624,23.072568284443342,40.05134539991026,22.324494864816362,59.93856433881809,50.509590303534864,32.84910148329592,55.89010710145792),c(39.758625325000324,40.9718859959571,40.19755261133803,51.46957540375822,39.1410875499123,34.41692958325856,22.30272161444394,46.901762424710945,60.03456477981668,39.92150507246994),c(83.13167113409159,79.53366105097555,75.3704111462588,94.06508608273055,120.15403619973078,78.13573202685727,101.75616736590047,70.95394828353707,98.54730444988776,103.79591318842185),c(21.68652290454563,41.77525630960332,24.40565694259809,37.27107184410078,49.15392389988986,31.626367725156513,52.968963834304354,49.30698101059356,49.8400160436214,55.89010710145792))
targetgene="DOT1L"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(54.21630726136408,44.98873756418819,63.16758267495976,84.30361488546606,62.80779164985927,48.36973887376879,69.69600504513731,43.293934545887026,58.901837142461645,49.50266628986273),c(68.67398919772783,77.1235501100369,105.51857560476233,74.54214368820156,93.75655854978993,116.27341075425188,68.30208494423457,113.0452735364828,126.86549538376354,92.61789176813026),c(22.590128025568365,24.10110940938653,13.638455350275404,33.72144595418642,33.67954044992454,34.41692958325856,32.06016232076316,28.86262303059135,24.9200080218107,31.937204057975954),c(35.24059971988665,24.904479723032747,28.712537579527165,46.145136568886684,30.038509049932696,37.2074914413606,54.3628839352071,40.88871596000442,37.380012032716046,51.09952649276153))
targetgene="NSD1"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(67.7703840767051,38.56177505501845,48.81131388519618,34.60885242666501,32.769282599926576,91.15835403133349,96.18048696228949,69.75133899059577,53.23819895568649,54.29324689855912),c(54.21630726136408,30.52807191855627,26.559097261062625,36.38366537162219,50.06418174988782,18.6037457206803,39.02976282527689,48.104371717652256,48.707288406266365,59.08382750725551),c(95.78214282840987,73.10669854180581,78.24166490421152,82.52880194050887,55.525728849875584,84.64704302909537,39.02976282527689,92.60091555648059,31.716373845940886,57.48696730435672),c(49.6982816562504,38.56177505501845,45.22224668775528,57.68142071110835,38.23082969991434,58.601799020142955,68.30208494423457,40.88871596000442,20.38909747239057,49.50266628986273))
targetgene="SP140"
collabel=c("day7_r1","day7_r2","day7_r3","day7_r4","day7_r5","day7_Cas9_r1","day7_Cas9_r2","day7_Cas9_r3","day7_Cas9_r4","day7_Cas9_r5")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("Day_7_summary.Rnw");
library(tools);

texi2dvi("Day_7_summary.tex",pdf=TRUE);

