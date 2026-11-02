pdf(file='Day7.5_Total.vs.Day7.5_Ly108_test.pdf',width=4.5,height=4.5);
gstable=read.table('Day7.5_Total.vs.Day7.5_Ly108_test.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Prdm1","Runx3","Ikzf1","Creb1","Gfi1","Rbpj","Nfat5","Irf2","Lhx2","Arnt")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Day7.5_Total_1,Day7.5_Total_2,Day7.5_Total_3_vs_Day7.5_Ly108_1,Day7.5_Ly108_2 neg.'


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
targetmat=list(c(621.0505063645882,719.5581822665099,50.736182372417986,82.3648078228476,86.08469283361052),c(1463.8147837639988,703.8086187695685,228.87655603557445,80.66656436258269,51.991745176735066),c(1444.9187237326219,1177.2798713963691,124.0217791325773,130.7647464403972,103.98349035347013),c(657.5828890919169,773.6973067872459,42.84388733670852,92.55426858443698,58.81033470811016))
targetgene="Prdm1"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(353.98619125446095,278.5704043521509,13.529648632644797,46.70169515728472,27.274358125500363),c(418.2327953611425,307.1164881903572,80.05042107648171,77.2700774420529,80.11842699365732),c(980.0756469607496,672.3094917756857,556.9705353772108,210.58218907284746,216.49021762115913),c(246.9085177433251,106.3095536043544,28.18676798467666,38.21047785596022,39.206889805406774))
targetgene="Runx3"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(3088.875946462414,1287.5268158749589,1769.001558718307,362.5749787665559,658.8462134691182),c(619.790769029163,724.4799208593041,296.52479919879846,220.77164983443686,238.65063359812817),c(233.0514070536487,377.98952392659345,92.45259898973944,168.12610256622497,115.06369834195466),c(28.973958714777943,88.59129467029534,14.657119352031863,5.943852110927146,2.5569710742656593))
targetgene="Ikzf1"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(211.6358723514215,188.99476196329672,180.3953151019306,133.31211163079456,121.02996418190786),c(444.68727940507017,395.7077828606525,241.27873394883218,173.2208329470197,182.39726996428368),c(207.8566603451461,184.07302337050254,151.0810763978669,83.21392955298005,104.83581404489202),c(500.11572216377584,249.03997279538578,145.44372280093157,174.91907640728456,173.8740330500648))
targetgene="Creb1"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(584.5181236372594,343.5373537770341,226.62161459680033,126.51913778973496,86.9370165250324),c(1970.2291926049002,899.693814762777,765.552618463818,331.1574747516553,138.92876170176748),c(114.63609752368664,96.46607641876604,60.883418846901584,90.85602512417209,103.13116666204824),c(18.89606003137692,7.874781748470697,19.167002229580127,11.038582491721842,10.227884297062637))
targetgene="Gfi1"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(534.1286302202543,400.62952145344667,189.41508085702714,146.89805931291374,181.5449462728618),c(624.8297183708635,582.7338493868316,294.2698577600243,193.59975447019846,191.77283056992442),c(617.2712943583127,613.2486286621555,463.3904656680843,534.9466899834431,518.2128043845069),c(23.93500937307743,0.0,31.56918014283786,24.624530173841034,24.717387051234706))
targetgene="Rbpj"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(3421.4466030146477,1243.2311685398113,1227.8156134125152,450.8836387003306,125.29158263901729),c(360.2848779315866,414.4103895132704,289.75997488247606,340.49781378311224,232.68436775817497),c(996.4522323212763,637.8573216261265,1518.7030590143784,425.40998679635715,330.7015922716919),c(32.75317072105333,74.81042661047162,91.32512827035238,73.02446879139065,79.26610330223544))
targetgene="Nfat5"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(180.1424389657933,214.58780264582649,3198.6344309011074,89.15778166390719,311.9504710604104),c(175.1034896240928,172.2608507477965,124.0217791325773,191.90151100993356,173.8740330500648),c(511.45335818260196,355.34952639974017,162.35578359173755,200.39272831125805,138.07643801034558),c(621.0505063645882,487.2521206866243,313.43685998960444,208.88394561258255,126.99623002186107))
targetgene="Irf2"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(4600.560748972567,1893.8850105072024,1952.779285978399,634.2939324089397,261.6633732665191),c(254.46694175587587,190.9634574004144,193.9249637345754,169.82434602648988,144.0427038502988),c(119.67504686538716,329.7564857172104,249.17102898454166,225.86638021523154,256.5494311179878),c(31.4934333856282,11.812172622706045,18.03953151019306,8.491217301324493,10.227884297062637))
targetgene="Lhx2"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(170.06454028239227,95.4817287002072,115.00201337748076,67.92973841059595,87.78934021645429),c(264.54484043927687,377.98952392659345,134.16901560706089,115.48055529801312,155.9752355302052),c(387.9990993109394,440.00343019580015,248.0435582651546,222.46989329470173,270.186610180738),c(192.7398123200446,30.51477927532395,107.1097183417713,44.15432996688737,69.03821900517279))
targetgene="Arnt"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetgenelist=c("Mnt","Fli1","Nfatc3","Hoxa10","Tcf7l2","Tgif1","Bcl6","Hltf","Zeb2","Tcf12")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Day7.5_Total_1,Day7.5_Total_2,Day7.5_Total_3_vs_Day7.5_Ly108_1,Day7.5_Ly108_2 pos.'


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
targetmat=list(c(280.92142579980356,87.6069469517365,183.7777272600918,208.88394561258255,164.49847244442407),c(325.012232539683,118.12172622706045,338.2412158161199,431.3538389072843,403.14910604255226),c(11.337636018826153,12.796520341264882,62.01088956628865,73.8735905215231,65.62892423948524),c(18.89606003137692,40.35825646091232,91.32512827035238,50.947303807946966,82.67539806792297))
targetgene="Mnt"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(773.4787239510285,544.3442883630369,821.9261544331714,695.430696978476,802.8889173194169),c(579.4791742955589,217.540845801503,564.8628304129203,444.93978658940347,544.6348388185854),c(148.6490055801651,241.16519104691508,629.128661417983,359.1784918460261,455.9931749107092),c(79.36345213178306,55.12347223929488,243.53367538760634,166.4278591059601,225.013454535378))
targetgene="Fli1"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(191.48007498461945,339.5999629027988,230.00402675496153,236.05584097682095,194.3298016441901),c(317.45380852713225,257.89910226241534,450.98828775482656,558.7220984271518,441.5036721565371),c(17.63632269595179,62.013906269206736,111.61960121931956,152.84191142384088,119.32531679906408),c(12.59737335425128,3.9373908742353483,38.334004459160255,24.624530173841034,33.240623965453565))
targetgene="Nfatc3"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(125.9737335425128,83.66955607750116,330.3489207804104,190.20326754966868,256.5494311179878),c(122.19452153623742,109.26259676003092,144.3162520815445,102.74372934602638,122.73461156475163),c(15.116848025101536,0.0,3.3824121581611992,1.6982434602648988,4.2616184571094315),c(15.116848025101536,44.29564733514767,121.76683769380317,151.143667963576,135.51946693607994))
targetgene="Hoxa10"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(212.89560968684663,171.27650302923766,208.58208308660727,90.00690339403964,126.99623002186107),c(307.3759098437312,241.16519104691508,298.77974063757256,363.42410049668837,328.99694488884813),c(76.8439774609328,69.88868801767744,828.6909787494938,58.58939937913901,428.71881678520884),c(148.6490055801651,23.62434524541209,21.42194366835426,73.8735905215231,41.76386087967243))
targetgene="Tcf7l2"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(118.41530952996203,72.84173117335395,196.17990517334954,225.86638021523154,205.41000963267462),c(5.038949341700512,25.593040682529764,72.15812604077225,45.852573427152265,64.77660054806336),c(129.75294554878818,105.32520588579557,170.24807862744703,163.03137218543029,141.48573277603313),c(27.714221379352814,41.34260417947116,43.97135805609559,113.78231183774822,74.15216115370411))
targetgene="Tgif1"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(318.7135458625574,189.97910968185556,312.30938927021737,599.4799414735093,461.1071170592405),c(61.72712943583127,111.23129219714859,152.20854711725394,180.01380678807928,134.66714324465804),c(12.59737335425128,36.420865586676975,23.676885107128392,15.284191142384088,22.160415976969045),c(23.93500937307743,30.51477927532395,49.608711653030916,116.32967702814557,89.49398759929807))
targetgene="Bcl6"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(132.27242021963843,143.71476690959022,127.4041912907385,163.03137218543029,190.92050687850255),c(327.5317072105333,326.8034425615339,197.30737589273662,279.3610492135759,225.86577822679988),c(51.64923075243025,22.639997526853254,54.11859453057919,50.09818207781451,40.05921349682866),c(11.337636018826153,2.9530431556765113,107.1097183417713,36.51223439569532,78.41377961081355))
targetgene="Hltf"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(141.09058156761432,114.1843353528251,154.4634885560281,213.1295542632448,205.41000963267462),c(221.7137710348225,182.10432793338487,127.4041912907385,166.4278591059601,141.48573277603313),c(11.337636018826153,46.264342772265344,146.57119352031862,61.985886299668806,91.19863498214184),c(18.89606003137692,5.906086311353023,1.1274707193870663,39.05959958609267,31.535976582609795))
targetgene="Zeb2"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
targetmat=list(c(200.29823633259537,244.1182342025916,347.2609815712164,273.4171971026487,321.32603166605116),c(142.35031890303947,139.77737603535488,147.6986642397057,172.37171121688723,204.5576859412527),c(102.03872416943537,57.09216767641255,85.68777467341704,135.8594768211919,156.8275592216271),c(277.1422137935282,282.50779522638624,588.5397155200486,413.52228257450287,495.20006471611595))
targetgene="Tcf12"
collabel=c("Day7.5_Ly108_1","Day7.5_Ly108_2","Day7.5_Total_1","Day7.5_Total_2","Day7.5_Total_3")

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
Sweave("Day7.5_Total.vs.Day7.5_Ly108_test_summary.Rnw");
library(tools);

texi2dvi("Day7.5_Total.vs.Day7.5_Ly108_test_summary.tex",pdf=TRUE);

