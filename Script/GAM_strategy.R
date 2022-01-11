library('mgcv')
library('mgcViz')
library(forcats)
library(itsadug)
library(gridExtra)
path=("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/")
file_name = paste(path,'GAM_table.csv',sep="")
data=read.csv(file_name,sep=';')                

#prepare and check data----
cols=c('Biome','Diadromous','Semelparity','Skip.breeding','Capital.income','Ovulation','Parental.care')
sapply(data[cols],class)
data$Biome=relevel(data$Biome,ref='Marine')
sapply(data[cols],levels)
levels(data$Semelparity)[levels(data$Semelparity)=="Iteroparous "] <-"Iteroparous"
levels(data$Diadromous)[levels(data$Diadromous)=="No "] <-"No"
sapply(data[cols],levels)
sapply(data,class)

#create ordered factors----
data$oBiome = as.ordered(data$Biome)
data$oDiadromous = as.ordered(data$Diadromous)
data$oSemelparity = as.ordered(data$Semelparity)
data$oSkip.breeding = as.ordered(data$Skip.breeding)
data$oCapital.income = as.ordered(data$Capital.income)
data$oOvulation = as.ordered(data$Ovulation)
data$oParental.care=as.ordered(data$Parental.care)
sapply(data,class)


## Automatic model selection with double penalty approach ##----
## Species as random intercept ##
## Repro.vars as ordered factors ##

m = bam(Log.λ.~ Biome + s(Environmental.autocorrelation,by=Biome)+
                Diadromous + s(Environmental.autocorrelation,by=Diadromous)+
                Semelparity + s(Environmental.autocorrelation,by=Semelparity)+
                Skip.breeding + s(Environmental.autocorrelation,by=Skip.breeding)+
                Capital.income + s(Environmental.autocorrelation,by=Capital.income)+
                Ovulation + s(Environmental.autocorrelation,by=Ovulation)+
                Parental.care + s(Environmental.autocorrelation,by=Parental.care)+
                Order + s(Environmental.autocorrelation,by=Order),
                data=data,
                select=TRUE,
                method="REML")
save(m, file="/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Reprostat_model.rda", compress='xz')
check_resid(m,split='Order')
summary(m)
(valRho <- acf(resid(m), plot=FALSE)$acf[2])
gam.check(m)

#With autocor structure
simdat <- start_event(data,column="Environmental.autocorrelation", event=c("Species"))
(valRho <- acf(resid(m), plot=FALSE)$acf[2])

m2 = bam(Log.λ.~ Biome + s(Environmental.autocorrelation,by=Biome)+
          Diadromous + s(Environmental.autocorrelation,by=Diadromous)+
          Semelparity + s(Environmental.autocorrelation,by=Semelparity)+
          Skip.breeding + s(Environmental.autocorrelation,by=Skip.breeding)+
          Capital.income + s(Environmental.autocorrelation,by=Capital.income)+
          Ovulation + s(Environmental.autocorrelation,by=Ovulation)+
          Parental.care + s(Environmental.autocorrelation,by=Parental.care)+
          Order + s(Environmental.autocorrelation,by=Order),
          data=simdat,
          AR.start=simdat$start.event, rho=valRho,
          select=TRUE,
          method="REML")
save(m2, file="/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Reprostat_model2.rda", compress='xz')
check_resid(m2,split='Order')
summary(m2)

#maintain relevant factors (Biome, Semelparity, Capital.income,Parental.care,Order)
# order relevant factors except 'Order'
m3 = bam(Log.λ.~s(Environmental.autocorrelation)+
                oBiome + s(Environmental.autocorrelation,by=oBiome)+
                oSemelparity + s(Environmental.autocorrelation,by=oSemelparity)+
                oSkip.breeding + s(Environmental.autocorrelation,by=oSkip.breeding)+
                oCapital.income + s(Environmental.autocorrelation,by=oCapital.income)+
                oParental.care + s(Environmental.autocorrelation,by=oParental.care)+
                Order + s(Environmental.autocorrelation,by=Order),
                data=simdat,
                AR.start=simdat$start.event,rho=valRho,
                select=TRUE,
                method="REML")
save(m3, file="/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Reprostat_model3.rda", compress='xz')
check_resid(m3)#,split='Order')
b=getViz(m3)
check(b)
summary(m3)
library(itsadug)
library(mgcv)
help(itsadug)
# select subset of
gamtabs(model=m3,type="HTML")
#explore plots
# Marine-vs-Freshwater
# Plot difference surface:
library('Cairo')

Cairo::Cairo(
  6, #length
  6, #width
  file = paste("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/test", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 1200,
  units = "cm" #you can change to pixels etc 
)


tiff("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Marine_Freshwater2.tiff", height = 20, width = 30, units='cm', 
     compression = "lzw", res = 300)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1.5, 0, 0)) 
p <- plot_diff(m3, view="Environmental.autocorrelation", 
               comp=list(oBiome=c("Marine", "Freshwater")),
               main="Marine-Freshwater",
               xlab =expression(bold('Environmental autocorrelation')),
               ylab=expression(bold("Difference smooth "*log*"("*lambda["s"]*")")),
               cex.lab=2.0,
               cex.axis=2.0,
               cex.main=2.5,
               col=rainbow(6))
dev.off()
# Iteroparous-vs-Semelparous
# Plot difference surface:
tiff("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Iteroparous_Semelparous2.tiff", height = 20, width = 30, units='cm', 
     compression = "lzw", res = 300)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1.5, 0, 0)) 
p2 <- plot_diff(m3, view="Environmental.autocorrelation", 
               comp=list(oSemelparity=c("Iteroparous", "Semelparous")),
               main="Iteroparous-Semelparous",
               xlab =expression(bold('Environmental autocorrelation')),
               ylab=expression(bold("Difference smooth "*log*"("*lambda["s"]*")")),
               cex.lab=2.0,
               cex.axis=2.0,
               cex.main=2.5,
               col=rainbow(6))
dev.off()

# Skip-vs-Obligate breeders
# Plot difference surface:
tiff("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Skip_breeding2.tiff", height = 20, width = 30, units='cm', 
     compression = "lzw", res = 300)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1.5, 0, 0)) 

p3 <- plot_diff(m3, view="Environmental.autocorrelation", 
                comp=list(oSkip.breeding=c("Yes", "No")),
                main="Skip - Obligate breeding",
                xlab =expression(bold('Environmental autocorrelation')),
                ylab=expression(bold("Difference smooth "*log*"("*lambda["s"]*")")),
                cex.lab=2.0,
                cex.axis=2.0,
                cex.main=2.5,
                col=rainbow(6))
dev.off()


# Capital-vs-Income
# Plot difference surface:
tiff("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Capital_Income2.tiff", height = 20, width = 30, units='cm', 
     compression = "lzw", res = 300)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1.5, 0, 0)) 

p4 <- plot_diff(m3, view="Environmental.autocorrelation", 
               comp=list(oCapital.income=c("Capital", "Income")),
               main="Capital - Income",
               xlab =expression(bold('Environmental autocorrelation')),
               ylab=expression(bold("Difference smooth "*log*"("*lambda["s"]*")")),
               cex.lab=2.0,
               cex.axis=2.0,
               cex.main=2.5,
               col=rainbow(6))
dev.off()


# Capital-vs-Income.capital
# Plot difference surface:
tiff("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Capital_income_capital2.tiff", height = 20, width = 30, units='cm', 
     compression = "lzw", res = 300)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1.5, 0, 0)) 

p5 <- plot_diff(m3, view="Environmental.autocorrelation", 
                comp=list(oCapital.income=c("Capital", "Income-Capital")),
                main="Capital - Income-Capital",
                xlab =expression(bold('Environmental autocorrelation')),
                ylab=expression(bold("Difference smooth "*log*"("*lambda["s"]*")")),
                cex.lab=2.0,
                cex.axis=2.0,
                cex.main=2.5,
                col=rainbow(6))
dev.off()

tiff("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/parental2.tiff", height = 20, width = 30, units='cm', 
     compression = "lzw", res = 300)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1.5, 0, 0)) 

# Parental-vs-Non-parental
# Plot difference surface:
p6 <- plot_diff(m3, view="Environmental.autocorrelation", 
                comp=list(oParental.care=c("Yes", "No")),
                main="Parental-investment - No Par. investment",
                xlab =expression(bold('Environmental autocorrelation')),
                ylab=expression(bold("Difference smooth "*log*"("*lambda["s"]*")")),
                cex.lab=2.0,
                cex.axis=2.0,
                cex.main=2.5,
                col=rainbow(6))
dev.off()

