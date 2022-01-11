library(ggplot2)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(wesanderson)
library(varrank)
library(lattice)
library(gridExtra) # also loads grid
library(latticeExtra)
library(car)
library(ggiraphExtra)
library(energy)


path=("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/")
file_name = paste(path,'varrank_table-Table 1.csv',sep="")
data=read.csv(file_name,sep=';')  

#drop the species column for varrank
data$Species <- NULL
df=as.data.frame(apply(apply(data, 2, gsub, patt=",", replace="."), 2, as.numeric))
sapply(df,class)

varrank.DEBIPMparam <- varrank(data.df =df, method = "estevez", variable.important = "Range.λ", discretization.method = "sturges", algorithm = "forward", scheme="mid", verbose = FALSE)
summary(varrank.DEBIPMparam)
plot(varrank.DEBIPMparam,main='MI Range λ | DEBIPM parameters',labelscex = 2,notecex = 2,maincex = 1.5)


# Plot relationships relevant variables to range(λ)
plat1 <- xyplot(df$Range.λ ~ log10(df$MuP), df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat2 <- xyplot(df$Range.λ ~ log10(df$Lp.cm), df)+layer(panel.smoother(..., col = "steelblue")) 
plat3 <- xyplot(df$Range.λ ~ log10(df$rB), df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat4 <- xyplot(df$Range.λ ~ log10(df$sigm.Lb.), df)+layer(panel.smoother(..., col = "steelblue")) 
plat5 <- xyplot(df$Range.λ ~ log10(df$Rm), df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat6 <- xyplot(df$Range.λ ~ log10(df$Phi), df)+layer(panel.smoother(..., col = "steelblue")) 
plat7 <- xyplot(df$Range.λ ~ log10(df$Lb..transformation.length.mm), df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat8 <- xyplot(df$Range.λ ~ log10(df$Lm.cm), df)+layer(panel.smoother(..., col = "steelblue")) 


#plat3 <- xyplot(df$Range.λ ~ df$rB, df)+layer(panel.smoother(..., col = "steelblue")) 

g=arrangeGrob(plat1, plat2,
             plat3,plat4,
             plat5,plat6,
             plat7,plat8,ncol = 3)               # Apply grid.arrange
g
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_log_smooths.png",dpi=300)

# Plot relationships relevant variables to range(λ)
plat1 <- xyplot(df$Range.λ ~ df$MuP, df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat2 <- xyplot(df$Range.λ ~ df$Lp.cm, df)+layer(panel.smoother(..., col = "steelblue")) 
plat3 <- xyplot(df$Range.λ ~ df$rB, df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat4 <- xyplot(df$Range.λ ~ df$sigm.Lb., df)+layer(panel.smoother(..., col = "steelblue")) 
plat5 <- xyplot(df$Range.λ ~ df$Rm, df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat6 <- xyplot(df$Range.λ ~ df$Phi, df)+layer(panel.smoother(..., col = "steelblue")) 
plat7 <- xyplot(df$Range.λ ~ df$Lb..transformation.length.mm, df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat8 <- xyplot(df$Range.λ ~ df$Lm.cm, df)+layer(panel.smoother(..., col = "steelblue")) 


#plat3 <- xyplot(df$Range.λ ~ df$rB, df)+layer(panel.smoother(..., col = "steelblue")) 

g=arrangeGrob(plat1, plat2,
              plat3,plat4,
              plat5,plat6,
              plat7,plat8,ncol = 3)               # Apply grid.arrange
g
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_untransformed_smooths.png",g,dpi=300)

#############
#REGRESSION#
#############
library(ggfortify)
df$Species = data$Species
df$logMuP=log10(df$MuP)
df$logLp=log10(df$Lp.cm)

#linear regresion-----
model =  lm(Range.λ ~ logMuP+logLp,data=df)
summary(model)
# diagnostic plots 
layout(matrix(c(1,2,3,4),2,2)) 
# optional 4 graphs/page 
autoplot(model)

# plot everything on one page
library(ggplot2)
library(jtools)
effect_plot(model, pred =logMuP,interval = TRUE,plot.points = TRUE)
effect_plot(model,pred=logLp,interval = TRUE,plot.points = TRUE)

#polynomial regression -----
modelpoly= lm(Range.λ  ~ poly(logMuP, 2) +logLp, data = df)
# diagnostic plots 
layout(matrix(c(1,2,3,4),2,2)) 
# optional 4 graphs/page 
autoplot(modelpoly)
effect_plot(modelpoly, pred =logMuP,interval = TRUE,plot.points = TRUE)
effect_plot(modelpoly,pred=logLp,interval = TRUE,plot.points = TRUE)
anova(model, modelpoly)



### drop outliers----
#4,17,21 (salmo salar, gillichthys mirabilis, anchoa mitchili)
clean_df = df[-c(17), ]

#scatterplot muP with species label/any clear clustering?
plot(Range.λ ~logMuP, col="lightblue", pch=19, cex=1,data=clean_df2)

### this add the labels to the points using values in speed

text(Range.λ ~logMuP, labels=Species,data=clean_df, cex=0.5, font=1)


#lm model
clean_model_lm = lm(Range.λ~logMuP+logLp,data=clean_df)
effect_plot(clean_model_lm, pred =logMuP,interval = TRUE,plot.points = TRUE)
effect_plot(clean_model_lm,pred=logLp,interval = TRUE,plot.points = TRUE)

######################################################################
### plot the model
######################################################################
layout(matrix(c(1,1),1,1)) 


summary(clean_model_lm)
anova(clean_model_lm,clean_model_lmm)
plot(clean_model_lm)
summary(model)
plat1 <- xyplot(clean_df$Range.λ ~ clean_df$logMuP, clean_df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat2 <- xyplot(clean_df$Range.λ ~ clean_df$logLp, clean_df)+layer(panel.smoother(..., col = "steelblue")) 
plat3 <- xyplot(clean_df$Range.λ ~ clean_df$MuP, clean_df)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat4 <- xyplot(clean_df$Range.λ ~ clean_df$Lp, clean_df)+layer(panel.smoother(..., col = "steelblue")) 




#grid.arrange(plat1,plat2,ncol=2)              # Apply grid.arrange

clean_df2 = df[-c(4,17,21), ]
clean_model_lm2 = lm(Range.λ~logMuP+logLp,data=clean_df2)
autoplot(clean_model_lm2)
summary(clean_model_lm2)

plat6 <- xyplot(clean_df2$Range.λ ~ clean_df2$logMuP, clean_df2)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat5 <- xyplot(clean_df2$Range.λ ~ clean_df2$logLp, clean_df2)+layer(panel.smoother(..., col = "steelblue")) 
plat7 <- xyplot(clean_df2$Range.λ ~ clean_df2$MuP, clean_df2)+layer(panel.smoother(..., col = "steelblue"))# Create plots
plat8 <- xyplot(clean_df2$Range.λ ~ clean_df2$Lp, clean_df2)+layer(panel.smoother(..., col = "steelblue")) 

grid.arrange(plat1,plat2,plat3,plat4,plat5,plat6,plat7,plat8,ncol=2)              # Apply grid.arrange

#fit polymodel
modelpoly2= lm(Range.λ  ~ poly(logMuP, 2) +poly(logLp,2), data = clean_df)
summary(modelpoly2)
effect_plot(modelpoly2, pred =logMuP,interval = TRUE,plot.points = TRUE)
effect_plot(modelpoly2,pred=logLp,interval = TRUE,plot.points = TRUE)



#                          Parabola y = a + bx + cx^2
#-------------------------------------------------------------------------------

model2 <- lm(Range.λ ~logMuP+I(logMuP^2)+logLp+I(logLp^2),data=clean_df2)
summary(model2)
coef(model2)
effect_plot(model2, pred =I(logMuP^2),interval = TRUE,plot.points = TRUE)
effect_plot(model2,pred=I(logLp^2),interval = TRUE,plot.points = TRUE)

# Predicted vs original
predicted <- fitted(model2)
original <- Range.λ


#                          Parabola y = a + bx + cx^2
#-------------------------------------------------------------------------------

model2 <- lm( ~logMuP+I(logMuP^2)+logLp+I(logLp^2),data=clean_df2)
model3= lm(Range.λ ~ poly(logMuP, 6)+poly(logLp,6), data =clean_df2) 
summary(model3)
model3

effect_plot(model3,pred=logMuP,interval=TRUE,plot.points=TRUE)
effect_plot(model3,pred=logLp,interval=TRUE,plot.points=TRUE)
# Predicted vs original
predicted <- fitted(model2)
original <- Range.λ

# Plot model2
# Plot model2
curve(predict(model2,clean_df2,col='green',lwd=2,add=TRUE))


#### range log lambda values shows 2 main groups <0.1(insensitive) and >0.1 (sensitive):
plot(cut(clean_df$Range.λ, breaks = 10),xlab='range λ',ylab='frequency')

##ggplot histogram of insensitivities 10 breaks
# Histogram for `chol$AGE`
hist= ggplot(data=clean_df, aes(Range.λ)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_histogram(binwidth=0.075,  
                # main = "Histogram sensitivity values", 
                 fill=I("gray"), 
                 col=I("black"), 
                 alpha=I(.5))+ 
  xlab("Sensitivity")+
  ylab("Frequency")+#+ylim(0,13)+
  scale_x_continuous(breaks = round(seq(0,1, by = 0.1),1))+
  scale_y_continuous(breaks = round(seq(0,12, by = 1),1))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=22,face="bold"))+
  geom_vline(xintercept = 0.115,linetype="dotted", color = "blue", size=2,alpha=0.75)
hist
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Sensitivity_histogram.png",hist,dpi=300)

#create new dataset with binary sensitivity
clean_df$bin_sens<-ifelse(clean_df$Range.λ<0.1,0,1)#,0,1)
clean_df$fac_sens<-ifelse(clean_df$Range.λ<0.1,'insensitive','sensitive')

#logistic regresiion
glm.fit <- glm(bin_sens ~ MuP + Lp.cm  + Phi + rB + sigm.Lb. +Rm +Lb..transformation.length.mm+Lm.cm, data = clean_df, family = binomial)
summary(glm.fit)

library(beanplot)
par(mfrow = c(4,2))
beanplot(MuP~bin_sens, ylab="MuP", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(Lp.cm~bin_sens, ylab="Lp", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(Phi~bin_sens, ylab="Phi", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(rB~bin_sens, ylab="rB", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(sigm.Lb.~bin_sens, ylab="Sigm(Lb)", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(Rm~bin_sens, ylab="Rm", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(Lb..transformation.length.mm~bin_sens, ylab="Lb", xlab= "Sens log λ", col="light blue",data = clean_df)
beanplot(Lm.cm~bin_sens, ylab="Lm", xlab= "Sens log λ", col="light blue",data = clean_df)



## conditional density plots
par(mfrow = c(4,2))
cdplot(as.factor(bin_sens) ~ log10(MuP), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(Lp.cm), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(Phi), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(rB), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(sigm.Lb.), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(Rm), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(Lb..transformation.length.mm), data=clean_df)
cdplot(as.factor(bin_sens) ~ log10(Lm.cm), data=clean_df)

table(as.factor(clean_df$bin_sens))  

## conditional density plots in ggplot----
library("ggsci")
thm <- theme_minimal() + theme(text = element_text(size = 16))

cp1= ggplot(clean_df) +
  geom_density(aes(x = log10(MuP), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(Natural mortality rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
cp1
cp2= ggplot(clean_df) +
  geom_density(aes(x = log10(Lp.cm), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(Maturation length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
cp3= ggplot(clean_df) +
  geom_density(aes(x = log10(Phi), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(Egg/larval survival)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
cp4= ggplot(clean_df) +
  geom_density(aes(x = log10(rB), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(von Bertallanfy growth rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

cp5= ggplot(clean_df) +
  geom_density(aes(x = log10(sigm.Lb.), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(variation in offspring size)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

cp6= ggplot(clean_df) +
  geom_density(aes(x = log10(Rm), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(Maximum reproduction rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
cp7= ggplot(clean_df) +
  geom_density(aes(x = log10(Lb..transformation.length.mm), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(Larval transformation length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
cp8= ggplot(clean_df) +
  geom_density(aes(x = log10(Lm.cm), y = ..count.., fill = fac_sens),alpha=0.7,
               position = "fill") +
  ylab('Conditional density') +
  xlab('Log10(Maximum length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
cpg=arrangeGrob(cp1, cp2,
              cp3,cp4,
              cp5,cp6,
              cp7,cp8,ncol = 2)               # Apply grid.arrange
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_conditional_density.png",cpg,dpi=300)

## density plots in ggplot----
library(wesanderson)
p1= ggplot(clean_df) +
  geom_density(aes(x = log10(MuP), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Natural mortality rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))


p2= ggplot(clean_df) +
  geom_density(aes(x = log10(Lp.cm), fill = fac_sens),alpha=0.7)+
  ylab('Density') +
  xlab('Log10(Maturation length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

p3= ggplot(clean_df) +
  geom_density(aes(x = log10(Phi), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Egg/larval survival)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

p4= ggplot(clean_df) +
  geom_density(aes(x = log10(rB), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(von Bertallanfy growth rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

p5= ggplot(clean_df) +
  geom_density(aes(x = log10(sigm.Lb.), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(variation in offspring size)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

p6= ggplot(clean_df) +
  geom_density(aes(x = log10(Rm), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Maximum reproduction rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

p7= ggplot(clean_df) +
  geom_density(aes(x = log10(Lb..transformation.length.mm), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Larval transformation length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

p8= ggplot(clean_df) +
  geom_density(aes(x = log10(Lm.cm), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Maximum length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))

pg=arrangeGrob(p1, p2,
                p3,p4,
                p5,p6,
                p7,p8,ncol = 2)               # Apply grid.arrange
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density.png",pg,dpi=300)




class(diamonds)
class(clean_df)

### Kolmogorov-smirnov test to see if parameter values between insensitive and sensitive species ----
####come from same distribution----

t1=ks.test(log10(clean_df$MuP[clean_df$fac_sens == "insensitive"]),log10(clean_df$MuP[clean_df$fac_sens == "sensitive"]))
t2=ks.test(log10(clean_df$Lp.cm[clean_df$fac_sens == "insensitive"]),log10(clean_df$Lp.cm[clean_df$fac_sens == "sensitive"]))
t3=ks.test(log10(clean_df$Phi[clean_df$fac_sens == "insensitive"]),log10(clean_df$Phi[clean_df$fac_sens == "sensitive"]))
t4=ks.test(log10(clean_df$rB[clean_df$fac_sens == "insensitive"]),log10(clean_df$rB[clean_df$fac_sens == "sensitive"]))
t5=ks.test(log10(clean_df$sigm.Lb.[clean_df$fac_sens == "insensitive"]),log10(clean_df$sigm.Lb.[clean_df$fac_sens == "sensitive"]))
t6=ks.test(log10(clean_df$Rm[clean_df$fac_sens == "insensitive"]),log10(clean_df$Rm[clean_df$fac_sens == "sensitive"]))
t7=ks.test(log10(clean_df$Lb..transformation.length.mm[clean_df$fac_sens == "insensitive"]),log10(clean_df$Lb..transformation.length.mm[clean_df$fac_sens == "sensitive"]))
t8=ks.test(log10(clean_df$Lm[clean_df$fac_sens == "insensitive"]),log10(clean_df$Lm[clean_df$fac_sens == "sensitive"]))

t1
t2
t3
t4
t5
t6
t7
t8
##### factorize DEB-Params for Fisher's exact test ###
clean_df$bin_MuP = cut(clean_df$MuP, breaks = c(0.33*max(clean_df$MuP),0.66*max(clean_df$MuP)))
plot(clean_df$bin_MuP)

library('Hmisc') # cut2
clean_df$bin_MuP <- as.factor(cut2(log10(clean_df$MuP), g=3))
clean_df$bin_Lp <- as.factor(cut2(log10(clean_df$Lp.cm), g=3))
clean_df$bin_Phi <- as.factor(cut2(log10(clean_df$Phi), g=3))
clean_df$bin_rB <- as.factor(cut2(log10(clean_df$rB), g=3))
clean_df$bin_sigmLb <- as.factor(cut2(log10(clean_df$sigm.Lb.), g=3))
clean_df$bin_Rm <- as.factor(cut2(log10(clean_df$Rm), g=3))
clean_df$bin_Lb <- as.factor(cut2(log10(clean_df$Lb..transformation.length.mm), g=3))
clean_df$bin_Lm <- as.factor(cut2(log10(clean_df$Lm.cm), g=3))

plot(clean_df$bin_rB)

library(dplyr)
library(timereg)
#names <- c("Low", "Medium", "High")
#b <- c(-Inf,(0.33*max(log10(clean_df$rB))), (0.66*max(log10(clean_df$rB))), Inf)
clean_df$bin_MuP = cut(log10(clean_df$MuP), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_Lp = cut(log10(clean_df$Lp.cm), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_Phi = cut(log10(clean_df$Phi), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_rB = cut(log10(clean_df$rB), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_sigmLb = cut(log10(clean_df$sigm.Lb.), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_Rm = cut(log10(clean_df$Rm), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_Lb = cut(log10(clean_df$Lb..transformation.length.mm), breaks=3)# ,labels=c('low','medium','high'))
clean_df$bin_Lm = cut(log10(clean_df$Lm.cm), breaks=3)# ,labels=c('low','medium','high'))

par(mfrow = c(4,2))
plot(clean_df$bin_MuP,main='MuP') 
plot(clean_df$bin_Lp,main='Lp') 
plot(clean_df$bin_Phi,main='Phi') 
plot(clean_df$bin_rB,main='rB') 
plot(clean_df$bin_sigmLb,main='Sigm(Lb)') 
plot(clean_df$bin_Rm,main='Rm') 
plot(clean_df$bin_Lb,main='Lb') 
plot(clean_df$bin_Lm,main='Lm')

library(rcompanion)
library(fmsb)

m1 <-  table(clean_df$fac_sens, clean_df$bin_MuP)
fisher.test(m1)
m2 <-  table(clean_df$fac_sens, clean_df$bin_Lp)
fisher.test(m2)
m3 <-  table(clean_df$fac_sens, clean_df$bin_Phi)
fisher.test(m3)
m4 <-  table(clean_df$fac_sens, clean_df$bin_rB)
fisher.test(m4)
m5 <-  table(clean_df$fac_sens, clean_df$bin_sigmLb)
fisher.test(m5)
m6 <-  table(clean_df$fac_sens, clean_df$bin_Rm)
fisher.test(m6)
m7 <-  table(clean_df$fac_sens, clean_df$bin_Lb)
fisher.test(m7)
m8 <-  table(clean_df$fac_sens, clean_df$bin_Lm)
fisher.test(m8)



## density plots in ggplot with fishers statistics----
library(wesanderson)
p1= ggplot(clean_df) +
  geom_density(aes(x = log10(MuP), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Natural mortality rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+
  geom_vline(xintercept = -0.542, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = -0.0847, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.169, p = 0.972", x = 0.2, y = 1.5,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 0.490", x = 0.2, y = 1.4,size=4,fontface='bold')
p1
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_1.png",p1,dpi=300)

p2= ggplot(clean_df) +
  geom_density(aes(x = log10(Lp.cm), fill = fac_sens),alpha=0.7)+
  ylab('Density') +
  xlab('Log10(Maturation length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = 0.85, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = 1.36, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
annotate("text",label="K-S test: D=0.287, p = 0.507", x = 1.65, y = 1.2,size=4,fontface='bold')+
annotate("text",label="Fisher's exact test: p = 0.240", x = 1.65, y = 1.1,size=4,fontface='bold')
p2
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_2.png",p2,dpi=300)

t3
p3= ggplot(clean_df) +
  geom_density(aes(x = log10(Phi), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Egg/larval survival)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = -2.06, linetype="dotted", 
           color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = -1.11, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.224, p = 0.801", x = -0.6, y = 0.5,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 1.000", x = -0.6, y = 0.45,size=4,fontface='bold')
p3
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_3.png",p3,dpi=300)

t4
p4= ggplot(clean_df) +
  geom_density(aes(x = log10(rB), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(von Bertallanfy growth rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = -0.546, linetype="dotted", 
           color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = -0.00109, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.290, p = 0.490", x = 0.3, y = 1.5,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 0.092", x = 0.3, y = 1.4,size=4,fontface='bold')
p4
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_4.png",p4,dpi=300)

t5
p5= ggplot(clean_df) +
  geom_density(aes(x = log10(sigm.Lb.), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(variation in offspring size)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = -0.324, linetype="dotted", 
           color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = 0.233, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.268, p = 0.593", x = 0.55, y = 0.9,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 0.469", x = 0.55, y = 0.83,size=4,fontface='bold')
p5
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_5.png",p5,dpi=300)


p6= ggplot(clean_df) +
  geom_density(aes(x = log10(Rm), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Maximum reproduction rate)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = 3.69, linetype="dotted", 
           color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = 5.23, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.290, p = 0.490", x = 6.1, y = 0.34,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 0.735", x = 6.1, y = 0.32,size=4,fontface='bold')

p6
t6
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_6.png",p6,dpi=300)


t7
p7= ggplot(clean_df) +
  geom_density(aes(x = log10(Lb..transformation.length.mm), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Larval transformation length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = 1.35, linetype="dotted", 
           color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = 1.63, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.232, p = 0.768", x = 1.8, y = 2.5,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 0.853", x = 1.8, y = 2.4,size=4,fontface='bold')
p7  
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_7.png",p7,dpi=300)


p8= ggplot(clean_df) +
  geom_density(aes(x = log10(Lm.cm), fill = fac_sens),alpha=0.7) +
  ylab('Density') +
  xlab('Log10(Maximum length)')+
  thm+ guides(fill=guide_legend(title="Response Log(λ)"))+scale_fill_manual(values=c("#999999","#00FF66"))+#scale_fill_manual(values=wes_palette(n = 2, "Royal1"))
  geom_vline(xintercept = 1.11, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  geom_vline(xintercept = 1.64, linetype="dotted", 
             color = "black", size=1,alpha=0.7)+
  annotate("text",label="K-S test: D=0.235, p = 0.751", x = 1.95, y = 0.85,size=4,fontface='bold')+
  annotate("text",label="Fisher's exact test: p = 0.295", x = 1.95, y = 0.8,size=4,fontface='bold')
p8
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test_8.png",p8,dpi=300)

t8

pg=arrangeGrob(p1, p2,
               p3,p4,
               p5,p6,
               p7,p8,ncol = 2)   # Apply grid.arrange
pg
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/varrank_table/DEB_param_density_test.png",pg,dpi=300)

##### Distance correlation for non-linear dependencies (all non-significant and close to zero, indicating independence)
dcor.test(log10(clean_df$MuP), clean_df$Range.λ,R=10000)#0.29
dcor.test(log10(clean_df$Lp.cm), clean_df$Range.λ,R=10000)
dcor.test(log10(clean_df$Phi), clean_df$Range.λ,R=10000)
dcor.test(log10(clean_df$rB), clean_df$Range.λ,R=10000)
dcor.test(log10(clean_df$sigm.Lb.), clean_df$Range.λ,R=10000)
dcor.test(log10(clean_df$Rm), clean_df$Range.λ,R=10000)
dcor.test(log10(clean_df$Lb..transformation.length.mm), clean_df$Range.λ,R=10000)
dcor.test(log10(clean_df$Lm.cm), clean_df$Range.λ,R=10000)

plot(clean_df$bin_MuP,main='MuP') 
plot(clean_df$bin_Lp,main='Lp') 
plot(clean_df$bin_Phi,main='Phi') 
plot(clean_df$bin_rB,main='rB') 
plot(clean_df$bin_sigmLb,main='Sigm(Lb)') 
plot(clean_df$bin_Rm,main='Rm') 
plot(clean_df$bin_Lb,main='Lb') 
plot(clean_df$bin_Lm,main='Lm')

#pearson correlation
cor.test(log10(clean_df$MuP), clean_df$Range.λ)#=10000)#0.29
cor.test(log10(clean_df$Lp.cm), clean_df$Range.λ)#,R=10000)
cor.test(log10(clean_df$Phi), clean_df$Range.λ)#,R=10000)
cor.test(log10(clean_df$rB), clean_df$Range.λ)#,R=10000)
cor.test(log10(clean_df$sigm.Lb.), clean_df$Range.λ)#,R=10000)
cor.test(log10(clean_df$Rm), clean_df$Range.λ)#,R=10000)
cor.test(log10(clean_df$Lb..transformation.length.mm), clean_df$Range.λ)#,R=10000)
cor.test(log10(clean_df$Lm.cm), clean_df$Range.λ)#,R=10000)




#manual distance correllation
doubleCenter <- function(x){
  centered <- x
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      centered[i,j] <- x[i,j] - mean(x[i,]) - mean(x[,j]) + mean(x)
    }
  }
  return(centered)
}

distanceCovariance <- function(x,y){
  N <- length(x)
  distX <- as.matrix(dist(x))
  distY <- as.matrix(dist(y))
  centeredX <- doubleCenter(distX)
  centeredY <- doubleCenter(distY)
  calc <- sum(centeredX * centeredY)
  return(sqrt(calc/(N^2)))
}

distanceVariance <- function(x){
  return(distanceCovariance(x,x))
}
distanceCorrelation <- function(x,y){
  cov <- distanceCovariance(x,y)
  sd <- sqrt(distanceVariance(x)*distanceVariance(y))
  return(cov/sd)
}

# Compare with Pearson's r
cor(log10(clean_df$MuP),clean_df$Range.λ) # 
distanceCorrelation(log10(clean_df$MuP),clean_df$Range.λ) # --> 0.509

bootstrap <- function(x,y,reps,alpha){
  estimates <- c()
  original <- data.frame(x,y)
  N <- dim(original)[1]
  for(i in 1:reps){
    S <- original[sample(1:N, N, replace = TRUE),]
    estimates <- append(estimates, distanceCorrelation(S$x, S$y))
  }
  u <- alpha/2 ; l <- 1-u
  interval <- quantile(estimates, c(l, u))
  return(2*(dcor(x,y)) - as.numeric(interval[1:2]))
}
bootstrap(log10(clean_df$MuP),clean_df$Range.λ,1000,0.05) 

permutationTest <- function(x,y,reps){
  estimates <- c()
  observed <- distanceCorrelation(x,y)
  N <- length(x)
  for(i in 1:reps){
    y_i <- sample(y, length(y), replace = T)
    estimates <- append(estimates, distanceCorrelation(x, y_i))
  }
  p_value <- mean(estimates >= observed)
  return(p_value)
}


permutationTest(log10(clean_df$MuP),clean_df$Range.λ,1000) # 
