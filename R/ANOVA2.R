apply_anova <- function(data) 
{
  print(deparse(substitute(data)))
    
  print("------------------------------------------------------")
  res.aov <- aov(as.formula(paste(colnames(data[1]), "~ Group * Condition")), data = data)
  
  print(summary(res.aov))
  
  # par(las=1,mar=c(3,5.5,0.1,1),tcl=(0.3), pty= "m", mgp=c(1.2,0.2,0),mfrow=c(1,1),cex.axis=0.8)
  # 
  # plot(TukeyHSD(res.aov,conf.level = 0.95))
  print("======================================================")
  return(res.aov)
}


library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(moments)
library(AID)
library(readxl)

meanData_Grip = read_excel("/Users/boramert/Desktop/Yüksek Lisans/Data_Results/Data/fNIRS_Data/Epoch_Means/meanData_Grip.xls");
meanData_Nback= read_excel("/Users/boramert/Desktop/Yüksek Lisans/Data_Results/Data/fNIRS_Data/Epoch_Means/meanData_Nback.xls");
meanData_Oddball = read_excel("/Users/boramert/Desktop/Yüksek Lisans/Data_Results/Data/fNIRS_Data/Epoch_Means/meanData_Oddball.xls");

for (roi in 1:12){
  data = cbind(meanData_Grip[roi],meanData_Grip[13],meanData_Grip[14])
  print(colnames(data[1]));
  res.aov <- apply_anova(data)
  print(summary(res.aov))
}

