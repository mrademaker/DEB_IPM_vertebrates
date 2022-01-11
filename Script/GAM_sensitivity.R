library(ggplot2)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(wesanderson)
#hrbrthemes::import_roboto_condensed()

path=("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/")
file_name = paste(path,'GAM_table.csv',sep="")
data=read.csv(file_name,sep=';')   
levels(data$Semelparity)[levels(data$Semelparity)=="Iteroparous "] <-"Iteroparous"

spec = unique(data$Species)
df_list=list()
df2=data.frame()
for (s in spec){
  print(s)
  spec_data = subset(data, Species == s)
  spec_data$range_Log.λ. = max(spec_data$Log.λ.) - min(spec_data$Log.λ.)
  df = data.frame(spec_data[!duplicated(spec_data[3:13]),])
  df2=rbind(df2,df)
  df_sens = df2[!duplicated(df2[1:13]),]
  df_sens$Environmental.autocorrelation = NULL
  df_sens$Log.λ. = NULL
  append(df_list,df_sens)
}

levels(df_sens$Semelparity)[levels(df_sens$Semelparity)=="Iteroparous "] <-"Iteroparous"
sapply(df_sens["Semelparity"],levels)

#Box plots and Tests for different factor levels
###BIOME### ----
# Plot 
df_sens_M_FW =  subset(df_sens, Biome == "Marine"|Biome=="Freshwater")

df_sens_M_FW%>%
   ggplot(aes(x=Biome, y=range_Log.λ., fill=Biome)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  geom_jitter(color="grey", size=7, alpha=0.6) +
  #theme_ipsum() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=32,face="bold"),
    axis.title.x=element_text(size=30,face="bold"),
    axis.title.y=element_text(size=30,face="bold"),
    axis.text.x=element_text(size=30),
    axis.text.y=element_text(size=30)) +
  ggtitle("Marine - Freshwater") +
  #scale_fill_manual(values=wes_palette(n = 2, "GrandBudapest1"))+
  scale_fill_manual(values=c("#666633", "#009999"))+
  xlab("Biome")+
  ylab(expression(Delta*log*"("*lambda["s"]*")"))

ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Marine_freshwater_sens2.png",dpi=300)
#Test
wilcox.test(df_sens_M_FW$range_Log.λ.~df_sens_M_FW$Biome,)
pairwise.wilcox.test(df_sens$range_Log.λ.,df_sens$Biome,p.adjust.method = "BH")#non sign.


###SEMELPARITY### ----
#Plot
df_sens %>%
  ggplot( aes(x=Semelparity, y=range_Log.λ., fill=Semelparity)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  geom_jitter(color="grey", size=7, alpha=0.6) +
  #theme_ipsum() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=32,face="bold"),
    axis.title.x=element_text(size=30,face="bold"),
    axis.title.y=element_text(size=30,face="bold"),
    axis.text.x=element_text(size=30),
    axis.text.y=element_text(size=30)) +
  #scale_fill_manual(values=c("99FF00", "#FF3333"))+
  scale_fill_manual(values=wes_palette(n = 2, "Chevalier1"))+
  ggtitle("Iteroparous - Semelparous") +
  xlab("Semelparity")+
  ylab(expression(Delta*log*"("*lambda["s"]*")"))
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Semelparity_sens2.png",dpi=300)

#Test
wilcox.test(df_sens$range_Log.λ.~df_sens$Semelparity)#sign.

#### Capital Income #### ----
#Plot
df_sens_C_I =  subset(df_sens, Capital.income == "Capital"|Capital.income=="Income"|Capital.income=="Income-Capital")
df_sens_C_I %>%
  ggplot( aes(x=Capital.income, y = range_Log.λ., fill=Capital.income)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  geom_jitter(color="grey", size=6.5, alpha=0.6) +
  #theme_ipsum() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=18,face="bold"),
    axis.title.x=element_text(size=14,face="bold"),
    axis.title.y=element_text(size=14,face="bold"),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14)) +
  ggtitle("Capital - Income - Income-Capital") +
  xlab("Breeding")+
  ylab(expression(Delta*log*"("*lambda["s"]*")"))

#Plot Capital - Income, Income-capital
df_sens_C_I =  subset(df_sens, Capital.income == "Capital"|Capital.income=="Income"|Capital.income=="Income-Capital")
df_sens_C_I %>%
  ggplot( aes(x=Capital.income, y = range_Log.λ., fill=Capital.income)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  geom_jitter(color="grey", size=5, alpha=0.6) +
  #theme_ipsum() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=18,face="bold"),
    axis.title.x=element_text(size=14,face="bold"),
    axis.title.y=element_text(size=14,face="bold"),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14)) +
  ggtitle("Capital - Income - Income-Capital") +
  xlab("Breeding")+
  ylab(expression(Delta*log*"("*lambda["s"]*")"))

#Test
kruskal.test(df_sens_C_I$range_Log.λ.~df_sens_C_I$Capital.income)
pairwise.wilcox.test(df_sens_C_I$range_Log.λ.,df_sens_C_I$Capital.income,p.adjust.method = "BH")#non sign.



####SKIP BREEDING #### ----
#Plot
df_sens %>%
  ggplot( aes(x=Skip.breeding, y = range_Log.λ., fill=Skip.breeding)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  geom_jitter(color="grey", size=7, alpha=0.6) +
  #theme_ipsum() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=32,face="bold"),
    axis.title.x=element_text(size=30,face="bold"),
    axis.title.y=element_text(size=30,face="bold"),
    axis.text.x=element_text(size=30),
    axis.text.y=element_text(size=30)) +
  scale_fill_manual(values=wes_palette(n = 2, "Moonrise3"))+
  ggtitle("Skip breeding") +
  xlab("Skip breeding")+
  ylab(expression(Delta*log*"("*lambda["s"]*")"))
ggsave("/Users/markrademaker/Documents/DEB-IPM/DEB_IPM_MATLAB_CODE/DataFrames_Posthoc/Skipbreeding_sens2.png",dpi=300)

#Test
wilcox.test(df_sens$range_Log.λ.~df_sens$Skip.breeding)#sign.

#### PARENTAL CARE #### ----
#Plot
df_sens %>%
  ggplot( aes(x=Parental.care, y = range_Log.λ., fill=Parental.care)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE) +
  geom_jitter(color="grey", size=5, alpha=0.6) +
  #theme_ipsum() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=18,face="bold"),
    axis.title.x=element_text(size=14,face="bold"),
    axis.title.y=element_text(size=14,face="bold"),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14)) +
  ggtitle("Parental care") +
  xlab("Parental care")+
  ylab(expression(Delta*log*"("*lambda["s"]*")"))

#Test
wilcox.test(df_sens$range_Log.λ.~df_sens$Parental.care)#non sign.
