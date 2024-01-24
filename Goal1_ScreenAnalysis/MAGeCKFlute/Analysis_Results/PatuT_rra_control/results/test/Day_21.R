pdf(file='Day_21.pdf',width=4.5,height=4.5);
gstable=read.table('Day_21.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("BRD1","EPC2","H2AFZ","KDM4E","KDM4B","JMJD8","TDRD10","TRIM33","SMARCD3","ARID4B")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='day21_Cas9_r1,day21_Cas9_r2,day21_Cas9_r3,day21_Cas9_r4,day21_Cas9_r5_vs_day21_r1,day21_r2,day21_r3,day21_r4,day21_r5 neg.'


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
targetmat=list(c(45.13221114344801,44.64464424528859,20.011594019155087,38.314382806367064,34.58664560071015,37.02076757646316,50.659130766024276,37.298163756006254,22.102056931493312,32.20117941918732),c(46.902101776524404,40.49165408293616,27.074509555327474,31.7729028150361,42.064839244106935,39.198459786843344,14.775579806757081,29.61913004153438,12.629746817996178,9.470935123290387),c(46.902101776524404,28.03268359587888,61.21193464682733,31.7729028150361,58.89077494174971,42.464998102413624,22.163369710135623,66.91729379754064,56.83386068098281,30.30699239452924),c(36.28275797806605,45.68289178587669,35.31457768086192,53.26633707226641,29.912774573587154,40.28730589203344,51.714529323649785,27.42512040882813,25.259493635992357,38.83083400549059))
targetgene="BRD1"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(92.03431291997242,62.2948524352864,75.3377657191721,95.31870844510831,111.23813044552723,63.15307410102538,88.65347884054249,94.34241420636877,95.77558003647103,92.81516420824579),c(4.424726582690981,10.382475405881067,2.354305178724128,11.213965699424506,10.282516259670585,1.0888461051900928,9.498587018629552,0.0,0.0,0.0),c(20.353742280378516,46.7211393264648,35.31457768086192,52.331839930647696,30.847548779011753,41.37615199722353,49.603732208398775,28.522125225181256,42.099156059987266,50.19595615343905),c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
targetgene="EPC2"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(114.15794583342732,146.39290322292302,69.45200277236178,114.00865127748249,78.52103325566628,75.1303812581164,88.65347884054249,129.44656832966876,88.40822772597325,71.0320134246779),c(46.902101776524404,50.874129488817225,29.4288147340516,42.98686851446061,64.4994201742973,25.043460419372135,79.15489182191294,47.171207103184386,63.1487340899809,58.7197977644004),c(47.7870470930626,25.956188514702667,36.49173027022398,30.838405673417395,18.69548410849197,50.08692083874427,25.329565383012138,32.91014449059376,18.944620226994267,28.412805369871162),c(23.008578229993105,49.835881948229115,58.8576294681032,45.79035993931674,53.28212970920212,30.487690945322598,77.04409470666192,26.328115592475005,27.364451438991722,45.46048859179386))
targetgene="H2AFZ"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(30.97308607883687,59.18010981352208,28.251662144689536,37.37988566474836,42.064839244106935,37.02076757646316,55.936123554151806,42.783187837771884,43.15163496148695,49.248862641110016),c(14.15912506461114,22.841445892938346,31.783119912775728,45.79035993931674,33.651871395285546,23.95461431418204,15.830978364382588,19.746086694356254,36.83676155248886,18.941870246580773),c(47.7870470930626,24.917940974114558,47.086103574482564,24.29692568208643,40.19529083325774,53.35345915431455,23.218768267761128,18.649081878003127,34.731803749489494,29.3598988822002),c(47.7870470930626,36.338663920583734,41.20034062767224,54.200834213885116,28.978000368162554,25.043460419372135,42.215942305020235,54.85024081765626,30.521888143490767,24.624431320555008))
targetgene="KDM4E"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(49.55693772613899,68.52433767881504,61.21193464682733,58.87331992197866,63.5646459688727,42.464998102413624,62.268514899904844,28.522125225181256,29.469409241991084,32.20117941918732),c(28.31825012922228,18.68845573058592,34.13742509149986,29.903908531798685,14.956387286793577,16.332691577851392,11.609384133880564,28.522125225181256,17.89214132549459,8.523841610961348),c(22.12363291345491,30.10917867705509,17.65728884043096,17.75544569075547,22.434580930190364,27.22115262975232,16.886376922008093,14.261062612590628,15.787183522495225,31.25408590685828),c(61.94617215767374,28.03268359587888,89.46359679151686,52.331839930647696,36.45619401155934,53.35345915431455,37.99434807451821,70.20830824660001,51.571466173484396,49.248862641110016))
targetgene="KDM4B"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(43.36232051037162,16.611960649409706,9.417220714896512,35.51089138151094,21.499806724765765,11.97730715709102,31.661956728765176,9.873043347178127,47.36155056748567,10.418028635619425),c(26.54835949614589,15.5737131088216,23.543051787241282,17.75544569075547,9.347742054245986,17.421537683041485,28.495761055888657,18.649081878003127,17.89214132549459,44.51339507946482),c(1.7698906330763926,2.0764950811762133,5.8857629468103205,0.0,0.9347742054245985,1.0888461051900928,2.1107971152510117,0.0,0.0,0.0),c(7.964507848843767,2.0764950811762133,11.771525893620641,9.34497141618709,0.0,1.0888461051900928,9.498587018629552,4.388019265412501,2.104957802999363,5.682561073974233))
targetgene="JMJD8"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(50.44188304267719,26.994436055290773,25.897356965965407,45.79035993931674,43.93438765495613,16.332691577851392,40.105145189769225,43.88019265412501,29.469409241991084,39.77792751781963),c(23.8935235465313,11.420722946469173,11.771525893620641,10.279468557805798,12.152064670519781,5.444230525950465,6.3323913457530345,3.2910144490593756,19.99709912849395,17.994776734251737),c(16.81396101422573,16.611960649409706,24.720204376603345,19.624439973992885,23.369355135614963,3.2665383155702785,5.276992788127529,15.358067428943752,11.577267916496497,14.206402684935581),c(46.01715645998621,29.070931136466985,12.948678482982704,25.23142282370514,39.26051662783314,50.08692083874427,16.886376922008093,20.843091510709378,26.31197253749204,14.206402684935581))
targetgene="TDRD10"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(59.29133620805915,49.835881948229115,67.09769759363765,52.331839930647696,42.99961344953153,45.7315364179839,45.382137977896754,44.97719747047813,43.15163496148695,60.61398478905848),c(57.52144557498276,31.1474262176432,57.68047687874114,30.838405673417395,64.4994201742973,31.57653705051269,52.769927881275294,28.522125225181256,36.83676155248886,31.25408590685828),c(36.28275797806605,29.070931136466985,32.96027250213779,24.29692568208643,31.78232298443635,26.132306524562228,29.551159613514162,37.298163756006254,37.889240453988535,16.10058970959366))
targetgene="TRIM33"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(73.4504612726703,77.868565544108,90.64074938087893,51.39734278902899,44.86916186038073,106.7069183086291,51.714529323649785,54.85024081765626,63.1487340899809,76.71457449865214),c(68.14078937344111,65.40959505705072,37.66888285958605,45.79035993931674,30.847548779011753,48.998074733554176,44.326739420271245,25.23111077612188,45.25659276448631,35.98955346850347),c(30.97308607883687,16.611960649409706,54.14901911065495,45.79035993931674,52.34735550377752,35.931921471273064,36.938949516892706,49.365216735890634,48.414029468985355,19.888963758909814),c(43.36232051037162,93.4422786529296,62.38908723618939,65.41479991330962,56.08645232547591,33.75422926089288,34.828152401641695,35.104154123300006,29.469409241991084,37.88374049316155))
targetgene="SMARCD3"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(103.53860203496897,90.32753603116528,92.99505455960306,60.742314205216076,85.99922689906306,96.90730336191827,101.31826153204855,95.43941902272189,112.61524246046594,115.54540850414273),c(29.20319544576048,55.02711965116965,43.554645806396366,46.72485708093544,38.32574242240854,28.309998734942415,48.548333650773266,26.328115592475005,41.04667715848758,19.888963758909814),c(50.44188304267719,79.94506062528421,68.27485018299971,106.53267414453282,90.67309792618606,43.55384420760372,67.54550768803237,55.947245634009384,38.94171935548822,63.45526532604559),c(38.93759392768064,36.338663920583734,34.13742509149986,35.51089138151094,51.41258129835292,51.17576694393436,36.938949516892706,52.65623118495001,53.67642397648376,51.14304966576809))
targetgene="ARID4B"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetgenelist=c("SUV39H2","PRKCD","BRWD1","ASXL3","CARM1","HCFC1","BAZ1B","CBX6","FBXO17","KMT2B")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='day21_Cas9_r1,day21_Cas9_r2,day21_Cas9_r3,day21_Cas9_r4,day21_Cas9_r5_vs_day21_r1,day21_r2,day21_r3,day21_r4,day21_r5 pos.'


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
targetmat=list(c(49.55693772613899,42.56814916411237,51.794713931930815,47.65935422255416,50.47780709292832,47.90922862836408,40.105145189769225,26.328115592475005,37.889240453988535,42.61920805480674),c(26.54835949614589,28.03268359587888,38.846035448948115,27.100417106942558,13.08683887594438,39.198459786843344,30.606558171139667,30.716134857887504,16.839662423994906,26.518618345213085),c(2.654835949614589,2.0764950811762133,9.417220714896512,0.0,0.0,0.0,3.1661956728765173,4.388019265412501,2.104957802999363,0.0),c(0.0,0.0,2.354305178724128,10.279468557805798,0.0,13.066153262281114,13.720181249131576,6.582028898118751,14.734704620995542,9.470935123290387))
targetgene="SUV39H2"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(20.353742280378516,19.726703271174028,12.948678482982704,14.017457124280634,6.54341943797219,6.533076631140557,18.997174037259104,16.45507224529688,29.469409241991084,18.941870246580773),c(46.01715645998621,59.18010981352208,44.73179839575843,60.742314205216076,35.521419806134745,58.797689680265016,50.659130766024276,65.82028898118752,75.77848090797707,60.61398478905848),c(15.929015697687534,23.87969343352645,23.543051787241282,18.68994283237418,20.56503251934117,40.28730589203344,29.551159613514162,26.328115592475005,23.154535832992995,28.412805369871162),c(46.902101776524404,31.1474262176432,34.13742509149986,36.445388523129644,31.78232298443635,27.22115262975232,33.77275384401619,46.07420228683126,63.1487340899809,48.301769128780975))
targetgene="PRKCD"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(25.663414179607692,23.87969343352645,34.13742509149986,33.64189709827352,28.978000368162554,51.17576694393436,22.163369710135623,63.62627934848126,28.416930340491405,51.14304966576809),c(33.62792202845146,33.22392129881941,44.73179839575843,57.93882278035995,33.651871395285546,54.44230525950464,34.828152401641695,46.07420228683126,57.886339582482485,59.66689127672944),c(61.94617215767374,33.22392129881941,47.086103574482564,47.65935422255416,26.17367775188876,77.30807346849659,37.99434807451821,31.81313967424063,32.62684594649013,24.624431320555008),c(74.33540658920849,96.55702127469392,76.51491830853416,123.35362269366958,124.32496932147161,89.2853806255876,117.14923989643114,121.7675346151969,111.56276355896625,86.18550962194253))
targetgene="BRWD1"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(48.6719924096008,59.18010981352208,50.617561342568756,57.93882278035995,45.80393606580533,43.55384420760372,46.437536535522256,26.328115592475005,23.154535832992995,45.46048859179386),c(69.02573468997932,66.44784259763883,89.46359679151686,62.611308488453496,88.80354951533687,106.7069183086291,116.09384133880565,111.89449126801877,101.03797454396943,73.87329396166501),c(46.902101776524404,43.60639670470048,45.9089509851205,29.903908531798685,56.08645232547591,20.688075998611765,30.606558171139667,35.104154123300006,19.99709912849395,38.83083400549059),c(32.742976711913265,13.497218027645387,37.66888285958605,14.951954265899342,15.891161492218176,30.487690945322598,35.8835509592672,29.61913004153438,12.629746817996178,29.3598988822002))
targetgene="ASXL3"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(51.326828359215384,45.68289178587669,48.26325616384462,22.42793139884901,39.26051662783314,51.17576694393436,68.60090624565788,75.69333232836564,26.31197253749204,45.46048859179386),c(71.6805706395939,73.71557538175557,85.93213902343068,76.62876561273413,49.54303288750372,109.97345662419937,50.659130766024276,96.53642383907501,55.78138177948313,91.86807069591676),c(51.326828359215384,37.37691146117184,40.023188038310174,32.70739995665481,37.39096821698394,56.61999746988483,53.825326438900795,39.49217338871251,36.83676155248886,32.20117941918732),c(7.07956253230557,10.382475405881067,18.834441429793024,8.41047427456838,30.847548779011753,20.688075998611765,44.326739420271245,29.61913004153438,49.46650837048504,28.412805369871162))
targetgene="CARM1"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(75.22035190574668,56.06536719175776,67.09769759363765,58.87331992197866,90.67309792618606,55.531151364694736,65.43471057278136,64.72328416483438,66.30617079447994,66.2965458630327),c(30.088140762298675,43.60639670470048,35.31457768086192,18.68994283237418,27.108451957313356,76.2192273633065,62.268514899904844,63.62627934848126,54.728902877983444,59.66689127672944),c(37.167703294604244,40.49165408293616,28.251662144689536,26.165919965323848,21.499806724765765,39.198459786843344,28.495761055888657,32.91014449059376,26.31197253749204,22.73024429589693),c(96.4590395026634,92.40403111234149,102.41227527449956,93.44971416187089,98.15129156958285,116.50653325533993,99.20746441679755,82.27536122648439,57.886339582482485,104.18028635619426))
targetgene="HCFC1"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(61.06122684113554,40.49165408293616,45.9089509851205,99.05669701158314,46.73871027122993,79.48576567887677,81.26568893716394,58.14125526671563,86.3032699229739,51.14304966576809),c(40.70748456075703,46.7211393264648,56.50332428937907,46.72485708093544,73.84716222854328,39.198459786843344,48.548333650773266,28.522125225181256,57.886339582482485,59.66689127672944),c(6.194617215767374,15.5737131088216,16.480136251068895,14.951954265899342,17.760709903067372,14.154999367471207,27.440362498263152,12.067052979884377,21.049578029993633,31.25408590685828))
targetgene="BAZ1B"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(78.76013317189947,25.956188514702667,31.783119912775728,69.15278847978446,26.17367775188876,45.7315364179839,62.268514899904844,59.23826008306876,65.25369189298026,55.87851722741328),c(35.39781266152785,56.06536719175776,44.73179839575843,28.034914248561268,23.369355135614963,69.68615073216594,17.9417754796336,28.522125225181256,45.25659276448631,43.56630156713578),c(54.86660962536817,48.79763440764101,64.74339241491352,64.48030277169092,52.34735550377752,75.1303812581164,59.102319227028325,58.14125526671563,38.94171935548822,73.87329396166501),c(21.238687596916712,17.650208189997812,27.074509555327474,0.9344971416187089,15.891161492218176,21.77692210380186,5.276992788127529,5.485024081765626,17.89214132549459,14.206402684935581))
targetgene="CBX6"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(84.95475038766685,87.21279340940096,78.8692234872583,58.87331992197866,67.30374279057109,64.24192020621548,59.102319227028325,68.01429861389376,73.67352310497772,62.50817181371656),c(3.539781266152785,0.0,5.8857629468103205,0.9344971416187089,0.0,8.710768841520743,0.0,3.2910144490593756,0.0,4.735467561645193),c(0.0,0.0,0.0,0.0,6.54341943797219,4.355384420760371,7.387789903378541,4.388019265412501,0.0,3.788374049316155),c(0.0,3.1147426217643197,1.177152589362064,0.0,0.0,3.2665383155702785,0.0,0.0,0.0,0.0))
targetgene="FBXO17"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
targetmat=list(c(52.211773675753584,53.988872110581546,29.4288147340516,40.183377089604484,42.064839244106935,40.28730589203344,42.215942305020235,74.59632751201251,61.04377628698153,67.24363937536175),c(4.424726582690981,20.764950811762134,22.365899197879216,16.82094854913676,17.760709903067372,17.421537683041485,40.105145189769225,26.328115592475005,26.31197253749204,15.15349619726462))
targetgene="KMT2B"
collabel=c("day21_r1","day21_r2","day21_r3","day21_r4","day21_r5","day21_Cas9_r1","day21_Cas9_r2","day21_Cas9_r3","day21_Cas9_r4","day21_Cas9_r5")

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
Sweave("Day_21_summary.Rnw");
library(tools);

texi2dvi("Day_21_summary.tex",pdf=TRUE);

