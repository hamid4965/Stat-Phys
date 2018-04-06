# In this script I will post process the results of the previous steps

rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)

## Box plot the smooth R2 for all the scales
R2_Smooth<-matrix(nrow = 9,ncol = 10)
leaf_smooth <- load("DryLeaf_N_Smooth.rds")
R2_Smooth[1,]<-r2_pls_plain
R2_Smooth[2,]<-r2_rf_plain
R2_Smooth[3,]<-r2_svm_plain
rm(list=setdiff(ls(), "R2_Smooth"))
shrub_smooth <- load("Canopy_N_Smooth_f.rds")
R2_Smooth[4,]<-r2_pls_plain
R2_Smooth[5,]<-r2_rf_plain
R2_Smooth[6,]<-r2_svm_plain
rm(list=setdiff(ls(), "R2_Smooth"))
plot_Smooth <- load("Plot_N_smooth_f.rds")
R2_Smooth[7,]<-r2_pls_plain
R2_Smooth[8,]<-r2_rf_plain
R2_Smooth[9,]<-r2_svm_plain
rm(list=setdiff(ls(), "R2_Smooth"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM","PLS","RF","SVM")
Scale <- c( "Leaf","Leaf","Leaf","Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2_Smooth)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Leaf","Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.9, y = 1.05, label = "(a)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))

## Box plot the Derv1 R2 for all the scales
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2_Derv1<-matrix(nrow = 9,ncol = 10)

leaf_derv1 <- load("DryLeaf_N_Derv.rds")
R2_Derv1[1,]<-r2_pls_plain
R2_Derv1[2,]<-r2_rf_plain
R2_Derv1[3,]<-r2_svm_plain
rm(list=setdiff(ls(), "R2_Derv1"))
shrub_derv1 <- load("Canopy_N_Derv1_f.rds")
R2_Derv1[4,]<-r2_pls_plain
R2_Derv1[5,]<-r2_rf_plain
R2_Derv1[6,]<-r2_svm_plain
rm(list=setdiff(ls(), "R2_Derv1"))
plot_derv1 <- load("Plot_N_Derv_f.rds")
R2_Derv1[7,]<-r2_pls_plain
R2_Derv1[8,]<-r2_rf_plain
R2_Derv1[9,]<-r2_svm_plain
rm(list=setdiff(ls(), "R2_Derv1"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM","PLS","RF","SVM")
Scale <- c( "Leaf","Leaf","Leaf","Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2_Derv1)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Leaf","Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.9, y = 1.05, label = "(b)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))

## Box plot the Log R2 for all the scales
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2_Log<-matrix(nrow = 9,ncol = 10)
leaf_log <- load("DryLeaf_N_Log.rds")
R2_Log[1,]<-r2_pls_plain
R2_Log[2,]<-r2_rf_plain
R2_Log[3,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2_Log"))
shrub_log <- load("Canopy_N_log_f.rds")
R2_Log[4,]<-r2_pls_plain
R2_Log[5,]<-r2_rf_plain
R2_Log[6,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2_Log"))
plot_log <- load("Plot_N_log.rds")
R2_Log[7,]<-r2_pls_plain
R2_Log[8,]<-r2_rf_plain
R2_Log[9,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2_Log"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM","PLS","RF","SVM")
Scale <- c( "Leaf","Leaf","Leaf","Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2_Log)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Leaf","Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.9, y = 1.05, label = "(c)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))

#------------------------------------------------------------------------------------------------------------------
## Box plot the Log R2 for all the scales
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2_dervLog<-matrix(nrow = 9,ncol = 10)
leaf_dervlog <- load("DryLeaf_N_DervLog.rds")
R2_dervLog[1,]<-r2_pls_plain
R2_dervLog[2,]<-r2_rf_plain
R2_dervLog[3,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2_dervLog"))
shrub_log <- load("Canopy_N_DervLog.rds")
R2_dervLog[4,]<-r2_pls_plain
R2_dervLog[5,]<-r2_rf_plain
R2_dervLog[6,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2_dervLog"))
plot_dervlog <- load("Plot_N_DervLog.rds")
R2_dervLog[7,]<-r2_pls_plain
R2_dervLog[8,]<-r2_rf_plain
R2_dervLog[9,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2_dervLog"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM","PLS","RF","SVM")
Scale <- c( "Leaf","Leaf","Leaf","Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2_dervLog)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Leaf","Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.9, y = 1.05, label = "(d)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))

#-------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2<-matrix(nrow = 6,ncol = 10)
canopy_smooth <- load("Corrected_Canopy_LAI_Smooth.rds")
R2[1,]<-r2_pls_plain
R2[2,]<-r2_rf_plain
R2[3,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))
plot_smooth <- load("Corrected_Plot_LAI_Smooth.rds")
R2[4,]<-r2_pls_plain
R2[5,]<-r2_rf_plain
R2[6,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM")
Scale <- c("Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.5, y = 1.05, label = "(a)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))



#-------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2<-matrix(nrow = 6,ncol = 10)
canopy_derv <- load("Corrected_Canopy_LAI_Derv1.rds")
R2[1,]<-r2_pls_plain
R2[2,]<-r2_rf_plain
R2[3,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))
plot_derv <- load("Corrected_Plot_LAI_Derv1.rds")
R2[4,]<-r2_pls_plain
R2[5,]<-r2_rf_plain
R2[6,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM")
Scale <- c("Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.5, y = 1.05, label = "(b)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))

#--------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2<-matrix(nrow = 6,ncol = 10)
canopy_log <- load("Corrected_Canopy_LAI_Log.rds")
R2[1,]<-r2_pls_plain
R2[2,]<-r2_rf_plain
R2[3,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))
plot_derv <- load("Corrected_Plot_LAI_Log.rds")
R2[4,]<-r2_pls_plain
R2[5,]<-r2_rf_plain
R2[6,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM")
Scale <- c("Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.5, y = 1.05, label = "(c)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))

#--------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()                            # close any open graphics
setwd("D:/StatPaper/final_working")              # set the directory 
#install.packages("gridExtra")
library("gridExtra")
library(ggplot2)
R2<-matrix(nrow = 6,ncol = 10)
canopy_dervlog <- load("Corrected_Canopy_LAI_Derv1Log.rds")
R2[1,]<-r2_pls_plain
R2[2,]<-r2_rf_plain
R2[3,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))
plot_dervlog <- load("Corrected_Plot_LAI_Derv1Log.rds")
R2[4,]<-r2_pls_plain
R2[5,]<-r2_rf_plain
R2[6,]<-r2_pls_plain
rm(list=setdiff(ls(), "R2"))

Method <- c("PLS","RF","SVM","PLS","RF","SVM")
Scale <- c("Canopy","Canopy","Canopy","Plot","Plot","Plot")
var1<-"R2_pls"
var2<-"R2_rf"
var3<-"R2_svm"
all.data<-data.frame(R2)
miny<-min(all.data)
maxy <- max(all.data)
all.data$Method <- Method
all.data$Scale <- Scale
stacked.data = melt(all.data, id = c("Scale","Method"))
stacked.data = stacked.data[, -3]
ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(1))+theme_light()+
  scale_x_discrete(name ="Scale",limits=c("Canopy","Plot"))+
  theme(axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title=element_text(size=14,face = "bold"),
        legend.position="top",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_grey(start = 0, end = 0.99)+labs(y="R2")+ ylim(0,1.05)+
  annotate("text",  x=1.5, y = 1.05, label = "(d)", size=7)+guides(fill=guide_legend(title="Ensemble fits:"))












