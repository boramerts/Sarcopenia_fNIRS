remove_outliers <- function(in_data) 
{
  # ggbetweenstats(data, Condition, theta, outlier.tagging = TRUE)
  Q <- quantile(in_data["theta"], probs=c(.25, .75), na.rm = TRUE)
  iqr <- IQR(unlist(in_data["theta"]))
  
  Lower <- Q[1] - 1.5*iqr
  Upper <- Q[2] + 1.5*iqr 
  
  eliminated<- subset(in_data, in_data["theta"] > (Q[1] - 1.5*iqr) & in_data["theta"] < (Q[2]+1.5*iqr))
  return(eliminated)
  # ggbetweenstats(eliminated, Condition, theta, outlier.tagging = TRUE)
}

apply_anova <- function(data) 
{
  print(deparse(substitute(data)))
  
  if(shapiro.test(data$theta)["p.value"]<0.05)
  {
    print("Data is not normal distributed!!!")
    sk = skewness(data$theta)
    print(paste(c("Skewness:",sk), collapse = " "))
    
    # if(sk<0)
    # {
    #   data$theta = sqrt(max(data$theta+1) - data$theta)
    # } else
    # {
    #   data$theta = sqrt(data$theta). 
    # } #BOX COX TRANSFORMATION YAP!!!
    
    print("------------------------------------------------------")
  }
  res.aov <- aov(rank(theta) ~ Group * Condition, data = data)
  
  print(summary(res.aov))
  
  # par(las=1,mar=c(3,5.5,0.1,1),tcl=(0.3), pty= "m", mgp=c(1.2,0.2,0),mfrow=c(1,1),cex.axis=0.8)
  # 
  # plot(TukeyHSD(res.aov,conf.level = 0.95))
  print("======================================================")
  return(res.aov)

}

Primary_Motor_Cortex_L_Grip <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Primary_Motor_Cortex_L_Nback <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Primary_Motor_Cortex_L_Oddball <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))

Primary_Motor_Cortex_R_Grip <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Primary_Motor_Cortex_R_Nback <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Primary_Motor_Cortex_R_Oddball <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))

Pre_Supplementary_Motor_Cortex_L_Grip <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Pre_Supplementary_Motor_Cortex_L_Nback <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Pre_Supplementary_Motor_Cortex_L_Oddball <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))

Pre_Supplementary_Motor_Cortex_R_Grip <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Pre_Supplementary_Motor_Cortex_R_Nback <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Pre_Supplementary_Motor_Cortex_R_Oddball <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))

Dorsolateral_Prefrontal_Cortex_L_Grip <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Dorsolateral_Prefrontal_Cortex_L_Nback <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Dorsolateral_Prefrontal_Cortex_L_Oddball <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))

Dorsolateral_Prefrontal_Cortex_R_Grip <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Dorsolateral_Prefrontal_Cortex_R_Nback <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))
Dorsolateral_Prefrontal_Cortex_R_Oddball <- matrix(0,ncol=3,dimnames = list(NULL,c("Group","Condition","theta")))

addToMatrix <- function(matrix,roi)
{
  if (grepl('Primary_Motor_Cortex_L', roi, fixed=TRUE)){
    if (grepl('Rest',matrix$Condition,fixed=TRUE)|grepl('Block',matrix$Condition,fixed=TRUE)){
       a = rbind(Primary_Motor_Cortex_L_Grip, matrix)
       Primary_Motor_Cortex_L_Grip <<- a
    } else if (grepl('Nback',matrix$Condition,fixed=TRUE)|grepl('Oback',matrix$Condition,fixed=TRUE)){
      a = rbind(Primary_Motor_Cortex_L_Nback, matrix)
      Primary_Motor_Cortex_L_Nback <<- a
    } else if (grepl('Odd',matrix$Condition,fixed=TRUE)|grepl('Std',matrix$Condition,fixed=TRUE)){
      a = rbind(Primary_Motor_Cortex_L_Oddball, matrix)
      Primary_Motor_Cortex_L_Oddball <<- a
    }
  } else if (grepl('Pre_Supplementary_Motor_Cortex_L', roi, fixed=TRUE)){
    if (grepl('Rest',matrix$Condition,fixed=TRUE)|grepl('Block',matrix$Condition,fixed=TRUE)){
      a = rbind(Pre_Supplementary_Motor_Cortex_L_Grip, matrix)
      Pre_Supplementary_Motor_Cortex_L_Grip <<- a
    } else if (grepl('Nback',matrix$Condition,fixed=TRUE)|grepl('Oback',matrix$Condition,fixed=TRUE)){
      a = rbind(Pre_Supplementary_Motor_Cortex_L_Nback, matrix)
      Pre_Supplementary_Motor_Cortex_L_Nback <<- a
    } else if (grepl('Odd',matrix$Condition,fixed=TRUE)|grepl('Std',matrix$Condition,fixed=TRUE)){
      a = rbind(Pre_Supplementary_Motor_Cortex_L_Oddball, matrix)
      Pre_Supplementary_Motor_Cortex_L_Oddball <<- a
    }
  } else if (grepl('Dorsolateral_Prefrontal_Cortex_L', roi, fixed=TRUE)){
    if (grepl('Rest',matrix$Condition,fixed=TRUE)|grepl('Block',matrix$Condition,fixed=TRUE)){
      a = rbind(Dorsolateral_Prefrontal_Cortex_L_Grip, matrix)
      Dorsolateral_Prefrontal_Cortex_L_Grip <<- a
    } else if (grepl('Nback',matrix$Condition,fixed=TRUE)|grepl('Oback',matrix$Condition,fixed=TRUE)){
      a = rbind(Dorsolateral_Prefrontal_Cortex_L_Nback, matrix)
      Dorsolateral_Prefrontal_Cortex_L_Nback <<- a
    } else if (grepl('Odd',matrix$Condition,fixed=TRUE)|grepl('Std',matrix$Condition,fixed=TRUE)){
      a = rbind(Dorsolateral_Prefrontal_Cortex_L_Oddball, matrix)
      Dorsolateral_Prefrontal_Cortex_L_Oddball <<- a
    }
  } else if (grepl('Primary_Motor_Cortex_R', roi, fixed=TRUE)){
    if (grepl('Rest',matrix$Condition,fixed=TRUE)|grepl('Block',matrix$Condition,fixed=TRUE)){
      a = rbind(Primary_Motor_Cortex_R_Grip, matrix)
      Primary_Motor_Cortex_R_Grip <<- a
    } else if (grepl('Nback',matrix$Condition,fixed=TRUE)|grepl('Oback',matrix$Condition,fixed=TRUE)){
      a = rbind(Primary_Motor_Cortex_R_Nback, matrix)
      Primary_Motor_Cortex_R_Nback <<- a
    } else if (grepl('Odd',matrix$Condition,fixed=TRUE)|grepl('Std',matrix$Condition,fixed=TRUE)){
      a = rbind(Primary_Motor_Cortex_R_Oddball, matrix)
      Primary_Motor_Cortex_R_Oddball <<- a
    }
  } else if (grepl('Pre_Supplementary_Motor_Cortex_R', roi, fixed=TRUE)){
    if (grepl('Rest',matrix$Condition,fixed=TRUE)|grepl('Block',matrix$Condition,fixed=TRUE)){
      a = rbind(Pre_Supplementary_Motor_Cortex_R_Grip, matrix)
      Pre_Supplementary_Motor_Cortex_R_Grip <<- a
    } else if (grepl('Nback',matrix$Condition,fixed=TRUE)|grepl('Oback',matrix$Condition,fixed=TRUE)){
      a = rbind(Pre_Supplementary_Motor_Cortex_R_Nback, matrix)
      Pre_Supplementary_Motor_Cortex_R_Nback <<- a
    } else if (grepl('Odd',matrix$Condition,fixed=TRUE)|grepl('Std',matrix$Condition,fixed=TRUE)){
      a = rbind(Pre_Supplementary_Motor_Cortex_R_Oddball, matrix)
      Pre_Supplementary_Motor_Cortex_R_Oddball <<- a
    }
  } else if (grepl('Dorsolateral_Prefrontal_Cortex_R', roi, fixed=TRUE)){
    if (grepl('Rest',matrix$Condition,fixed=TRUE)|grepl('Block',matrix$Condition,fixed=TRUE)){
      a = rbind(Dorsolateral_Prefrontal_Cortex_R_Grip, matrix)
      Dorsolateral_Prefrontal_Cortex_R_Grip <<- a
    } else if (grepl('Nback',matrix$Condition,fixed=TRUE)|grepl('Oback',matrix$Condition,fixed=TRUE)){
      a = rbind(Dorsolateral_Prefrontal_Cortex_R_Nback, matrix)
      Dorsolateral_Prefrontal_Cortex_R_Nback <<- a
    } else if (grepl('Odd',matrix$Condition,fixed=TRUE)|grepl('Std',matrix$Condition,fixed=TRUE)){
      a = rbind(Dorsolateral_Prefrontal_Cortex_R_Oddball, matrix)
      Dorsolateral_Prefrontal_Cortex_R_Oddball <<- a
    }
  }
}



sarco = c('006', '009', '014', '016', '019', '020', '023', '025','037', '038', '039')

filenames = list.files(path="/Users/boramert/Desktop/Yüksek Lisans/Exports/GLM_Data",pattern="*.csv", full.names=TRUE)

for (i in 1:length(filenames))
{
  is_sarco = FALSE
  data = read.csv(filenames[i])
  data = data[data$Chroma=='hbo',1:4]
  
  for (k in 1:length(sarco))
  {
    if (grepl(sarco[k], filenames[i], fixed=FALSE)){
      is_sarco = TRUE
      break
    }
  }
  
  if (!is_sarco)
  {
    for (k in 1:nrow(data))
    {
      roi = data[k,"ROI"]
      m = data[k,c(2,4)]
      m <- cbind(Group='Control',m)
      addToMatrix(m,roi)
    }
  } else
  {
    for (k in 1:nrow(data))
    {
      roi = data[k,"ROI"]
      m = data[k,c(2,4)]
      m <- cbind(Group='Sarco',m)
      # print("____________________________________")
      # print(roi)
      # print(m)
      # print("----------------------------------")
      addToMatrix(m,roi)
    }
  }
}

Primary_Motor_Cortex_L_Grip <- Primary_Motor_Cortex_L_Grip[-1,]
Primary_Motor_Cortex_L_Nback <- Primary_Motor_Cortex_L_Nback[-1,]
Primary_Motor_Cortex_L_Oddball <- Primary_Motor_Cortex_L_Oddball[-1,]

Primary_Motor_Cortex_R_Grip <- Primary_Motor_Cortex_R_Grip[-1,]
Primary_Motor_Cortex_R_Nback <- Primary_Motor_Cortex_R_Nback[-1,]
Primary_Motor_Cortex_R_Oddball <- Primary_Motor_Cortex_R_Oddball[-1,]

Pre_Supplementary_Motor_Cortex_L_Grip <- Pre_Supplementary_Motor_Cortex_L_Grip[-1,]
Pre_Supplementary_Motor_Cortex_L_Nback <- Pre_Supplementary_Motor_Cortex_L_Nback[-1,]
Pre_Supplementary_Motor_Cortex_L_Oddball <- Pre_Supplementary_Motor_Cortex_L_Oddball[-1,]

Pre_Supplementary_Motor_Cortex_R_Grip <- Pre_Supplementary_Motor_Cortex_R_Grip[-1,]
Pre_Supplementary_Motor_Cortex_R_Nback <- Pre_Supplementary_Motor_Cortex_R_Nback[-1,]
Pre_Supplementary_Motor_Cortex_R_Oddball <- Pre_Supplementary_Motor_Cortex_R_Oddball[-1,]

Dorsolateral_Prefrontal_Cortex_L_Grip <- Dorsolateral_Prefrontal_Cortex_L_Grip[-1,]
Dorsolateral_Prefrontal_Cortex_L_Nback <- Dorsolateral_Prefrontal_Cortex_L_Nback[-1,]
Dorsolateral_Prefrontal_Cortex_L_Oddball <- Dorsolateral_Prefrontal_Cortex_L_Oddball[-1,]

Dorsolateral_Prefrontal_Cortex_R_Grip <- Dorsolateral_Prefrontal_Cortex_R_Grip[-1,]
Dorsolateral_Prefrontal_Cortex_R_Nback <- Dorsolateral_Prefrontal_Cortex_R_Nback[-1,]
Dorsolateral_Prefrontal_Cortex_R_Oddball <- Dorsolateral_Prefrontal_Cortex_R_Oddball[-1,]

library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(moments)
library(AID)

# Remove Outliers
Primary_Motor_Cortex_L_Grip <- remove_outliers(Primary_Motor_Cortex_L_Grip)
Primary_Motor_Cortex_L_Nback <- remove_outliers(Primary_Motor_Cortex_L_Nback)
Primary_Motor_Cortex_L_Oddball <- remove_outliers(Primary_Motor_Cortex_L_Oddball)

Primary_Motor_Cortex_R_Grip <- remove_outliers(Primary_Motor_Cortex_R_Grip)
Primary_Motor_Cortex_R_Nback <- remove_outliers(Primary_Motor_Cortex_R_Nback)
Primary_Motor_Cortex_R_Oddball <- remove_outliers(Primary_Motor_Cortex_R_Oddball)

Pre_Supplementary_Motor_Cortex_L_Grip <- remove_outliers(Pre_Supplementary_Motor_Cortex_L_Grip)
Pre_Supplementary_Motor_Cortex_L_Nback <- remove_outliers(Pre_Supplementary_Motor_Cortex_L_Nback)
Pre_Supplementary_Motor_Cortex_L_Oddball <- remove_outliers(Pre_Supplementary_Motor_Cortex_L_Oddball)

Pre_Supplementary_Motor_Cortex_R_Grip <- remove_outliers(Pre_Supplementary_Motor_Cortex_R_Grip)
Pre_Supplementary_Motor_Cortex_R_Nback <- remove_outliers(Pre_Supplementary_Motor_Cortex_R_Nback)
Pre_Supplementary_Motor_Cortex_R_Oddball <- remove_outliers(Pre_Supplementary_Motor_Cortex_R_Oddball)

Dorsolateral_Prefrontal_Cortex_L_Grip <- remove_outliers(Dorsolateral_Prefrontal_Cortex_L_Grip)
Dorsolateral_Prefrontal_Cortex_L_Nback <- remove_outliers(Dorsolateral_Prefrontal_Cortex_L_Nback)
Dorsolateral_Prefrontal_Cortex_L_Oddball <- remove_outliers(Dorsolateral_Prefrontal_Cortex_L_Oddball)

Dorsolateral_Prefrontal_Cortex_R_Grip <- remove_outliers(Dorsolateral_Prefrontal_Cortex_R_Grip)
Dorsolateral_Prefrontal_Cortex_R_Nback <- remove_outliers(Dorsolateral_Prefrontal_Cortex_R_Nback)
Dorsolateral_Prefrontal_Cortex_R_Oddball <- remove_outliers(Dorsolateral_Prefrontal_Cortex_R_Oddball)

# Visualize Data
plot_1 <- ggboxplot(Primary_Motor_Cortex_L_Grip, x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_L_Grip)), ylim = c(-1e-7,2.5e-7))
plot_2 <- ggboxplot(Primary_Motor_Cortex_R_Grip, x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_R_Grip)), ylim = c(-1e-7,2.5e-7))
plot_3 <- ggboxplot(Pre_Supplementary_Motor_Cortex_L_Grip, x="Condition",y="theta",add="point", title = deparse(substitute(Pre_Supplementary_Motor_Cortex_L_Grip)), ylim = c(-1e-7,2.5e-7))
plot_4 <- ggboxplot(Pre_Supplementary_Motor_Cortex_R_Grip, x="Condition",y="theta",add="point", title = deparse(substitute(Pre_Supplementary_Motor_Cortex_R_Grip)), ylim = c(-1e-7,2.5e-7))
plot_5 <- ggboxplot(Dorsolateral_Prefrontal_Cortex_L_Grip, x="Condition",y="theta",add="point", title = deparse(substitute(Dorsolateral_Prefrontal_Cortex_L_Grip)), ylim = c(-1e-7,2.5e-7))
plot_6 <- ggboxplot(Dorsolateral_Prefrontal_Cortex_R_Grip, x="Condition",y="theta",add="point", title = deparse(substitute(Dorsolateral_Prefrontal_Cortex_R_Grip)), ylim = c(-1e-7,2.5e-7))
grip_plot <- ggarrange(plot_1,
                          plot_2, plot_3, plot_4, plot_5, plot_6,
                          nrow = 3,
                          ncol = 2)

plot_1 <- ggboxplot(Primary_Motor_Cortex_L_Nback, x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_L_Nback)), ylim = c(-1e-7,2.5e-7))
plot_2 <- ggboxplot(Primary_Motor_Cortex_R_Nback, x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_R_Nback)), ylim = c(-1e-7,2.5e-7))
plot_3 <- ggboxplot(Pre_Supplementary_Motor_Cortex_L_Nback, x="Condition",y="theta",add="point", title = deparse(substitute(Pre_Supplementary_Motor_Cortex_L_Nback)), ylim = c(-1e-7,2.5e-7))
plot_4 <- ggboxplot(Pre_Supplementary_Motor_Cortex_R_Nback, x="Condition",y="theta",add="point", title = deparse(substitute(Pre_Supplementary_Motor_Cortex_R_Nback)), ylim = c(-1e-7,2.5e-7))
plot_5 <- ggboxplot(Dorsolateral_Prefrontal_Cortex_L_Nback, x="Condition",y="theta",add="point", title = deparse(substitute(Dorsolateral_Prefrontal_Cortex_L_Nback)), ylim = c(-1e-7,2.5e-7))
plot_6 <- ggboxplot(Dorsolateral_Prefrontal_Cortex_R_Nback, x="Condition",y="theta",add="point", title = deparse(substitute(Dorsolateral_Prefrontal_Cortex_R_Nback)), ylim = c(-1e-7,2.5e-7))
nback_plot <- ggarrange(plot_1,
                           plot_2, plot_3, plot_4, plot_5, plot_6,
                           nrow = 3,
                           ncol = 2)

plot_1 <- ggboxplot(Primary_Motor_Cortex_L_Oddball, x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_L_Oddball)), ylim = c(-1e-7,2.5e-7))
plot_2 <- ggboxplot(Primary_Motor_Cortex_R_Oddball, x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_R_Oddball)), ylim = c(-1e-7,2.5e-7))
plot_3 <- ggboxplot(Pre_Supplementary_Motor_Cortex_L_Oddball, x="Condition",y="theta",add="point", title = deparse(substitute(Pre_Supplementary_Motor_Cortex_L_Oddball)), ylim = c(-1e-7,2.5e-7))
plot_4 <- ggboxplot(Pre_Supplementary_Motor_Cortex_R_Oddball, x="Condition",y="theta",add="point", title = deparse(substitute(Pre_Supplementary_Motor_Cortex_R_Oddball)), ylim = c(-1e-7,2.5e-7))
plot_5 <- ggboxplot(Dorsolateral_Prefrontal_Cortex_L_Oddball, x="Condition",y="theta",add="point", title = deparse(substitute(Dorsolateral_Prefrontal_Cortex_L_Oddball)), ylim = c(-1e-7,2.5e-7))
plot_6 <- ggboxplot(Dorsolateral_Prefrontal_Cortex_R_Oddball, x="Condition",y="theta",add="point", title = deparse(substitute(Dorsolateral_Prefrontal_Cortex_R_Oddball)), ylim = c(-1e-7,2.5e-7))
odd_plot <- ggarrange(plot_1,
                           plot_2, plot_3, plot_4, plot_5, plot_6,
                           nrow = 3,
                           ncol = 2)

# ANOVA
ANOVA_Primary_Motor_Cortex_L_Grip <- apply_anova(Primary_Motor_Cortex_L_Grip)
ANOVA_Primary_Motor_Cortex_R_Grip <- apply_anova(Primary_Motor_Cortex_R_Grip)
ANOVA_Pre_Supplementary_Motor_Cortex_L_Grip <- apply_anova(Pre_Supplementary_Motor_Cortex_L_Grip)
ANOVA_Pre_Supplementary_Motor_Cortex_R_Grip <- apply_anova(Pre_Supplementary_Motor_Cortex_R_Grip)
ANOVA_Dorsolateral_Prefrontal_Cortex_L_Grip <- apply_anova(Dorsolateral_Prefrontal_Cortex_L_Grip)
ANOVA_Dorsolateral_Prefrontal_Cortex_R_Grip <- apply_anova(Dorsolateral_Prefrontal_Cortex_R_Grip)

ANOVA_Primary_Motor_Cortex_L_Nback <- apply_anova(Primary_Motor_Cortex_L_Nback)
ANOVA_Primary_Motor_Cortex_R_Nback <- apply_anova(Primary_Motor_Cortex_R_Nback)
ANOVA_Pre_Supplementary_Motor_Cortex_L_Nback <- apply_anova(Pre_Supplementary_Motor_Cortex_L_Nback)
ANOVA_Pre_Supplementary_Motor_Cortex_R_Nback <- apply_anova(Pre_Supplementary_Motor_Cortex_R_Nback)
ANOVA_Dorsolateral_Prefrontal_Cortex_L_Nback <- apply_anova(Dorsolateral_Prefrontal_Cortex_L_Nback)
ANOVA_Dorsolateral_Prefrontal_Cortex_R_Nback <- apply_anova(Dorsolateral_Prefrontal_Cortex_R_Nback)

ANOVA_Primary_Motor_Cortex_L_Oddball <- apply_anova(Primary_Motor_Cortex_L_Oddball)
ANOVA_Primary_Motor_Cortex_R_Oddball <- apply_anova(Primary_Motor_Cortex_R_Oddball)
ANOVA_Pre_Supplementary_Motor_Cortex_L_Oddball <- apply_anova(Pre_Supplementary_Motor_Cortex_L_Oddball)
ANOVA_Pre_Supplementary_Motor_Cortex_R_Oddball <- apply_anova(Pre_Supplementary_Motor_Cortex_R_Oddball)
ANOVA_Dorsolateral_Prefrontal_Cortex_L_Oddball <- apply_anova(Dorsolateral_Prefrontal_Cortex_L_Oddball)
ANOVA_Dorsolateral_Prefrontal_Cortex_R_Oddball <- apply_anova(Dorsolateral_Prefrontal_Cortex_R_Oddball)

#Visualization
# ggboxplot(Primary_Motor_Cortex_L_Grip, x="Condition",y="theta",add="point", title = "Primary_Motor_Cortex_L_Grip")
# ggboxplot(Primary_Motor_Cortex_R_Grip, x="Condition",y="theta",add="point", title = "Primary_Motor_Cortex_R_Grip")

plot_1 <- ggboxplot(Primary_Motor_Cortex_L_Grip[Primary_Motor_Cortex_L_Grip$Group == 'Sarco',], x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_L_Grip)), ylim = c(-1e-7,2.5e-7))
plot_2 <- ggboxplot(Primary_Motor_Cortex_L_Grip[Primary_Motor_Cortex_L_Grip$Group == 'Control',], x="Condition",y="theta",add="point", title = deparse(substitute(Primary_Motor_Cortex_L_Grip)), ylim = c(-1e-7,2.5e-7))
a <- ggarrange(plot_1,
                      plot_2,
                      nrow = 1,
                      ncol = 2)