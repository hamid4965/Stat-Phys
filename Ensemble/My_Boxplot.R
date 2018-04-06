My_Boxplot <- function(){
  
  leaf_smooth <- readRDS("leaf.out.smooth.ensemble.rds")
  shrub_smooth <- readRDS("shrub.out.smooth.ensemble.rds")
  plot_smooth <- readRDS("plot.out.smooth.ensemble.rds")
  leaf_derv <- readRDS("leaf.out.derv.ensemble.rds")
  shrub_derv <- readRDS("shrub.out.derv.ensemble.rds")
  plot_derv <- readRDS("plot.out.derv.ensemble.rds")
  leaf_log <- readRDS("leaf.out.log.ensemble.rds")
  shrub_log <- readRDS("shrub.out.log.ensemble.rds")
  plot_log <- readRDS("plot.out.log.ensemble.rds")
  
  
  
  Method <- c("PLSR","RF","SVM","PLSR","RF","SVM","PLSR","RF","SVM")
  Scale <- c( "Leaf","Leaf","Leaf","Canopy","Canopy","Canopy","Plot","Plot","Plot")
  var1<-"RMSEPs_pls"
  var2<-"RMSEPs_rf"
  var3<-"RMSEPs_svm"
  
  # all.data<-data.frame(rbind(leaf_smooth[[var1]],leaf_smooth[[var2]],leaf_smooth[[var3]], # Smoothed dataset
  #                           shrub_smooth[[var1]],shrub_smooth[[var2]],shrub_smooth[[var3]],
  #                          plot_smooth[[var1]],plot_smooth[[var2]],plot_smooth[[var3]]))
  
  # all.data<-data.frame(rbind(leaf_derv[[var1]],leaf_derv[[var2]],leaf_derv[[var3]], # first derivative dataset
  #                            shrub_derv[[var1]],shrub_derv[[var2]],shrub_derv[[var3]],
  #                            plot_derv[[var1]],plot_derv[[var2]],plot_derv[[var3]]))
  all.data<-data.frame(rbind(leaf_log[[var1]],leaf_log[[var2]],leaf_log[[var3]], # first derivative dataset
                             shrub_log[[var1]],shrub_log[[var2]],shrub_log[[var3]],
                             plot_log[[var1]],plot_log[[var2]],plot_log[[var3]]))
  
  
  miny<-min(all.data)
  maxy <- max(all.data)
  all.data$Method <- Method
  all.data$Scale <- Scale
  
  stacked.data = melt(all.data, id = c("Scale","Method"))
  stacked.data = stacked.data[, -3]
  
  
  
  ggplot(stacked.data, aes(Scale, value)) + geom_boxplot(aes(fill = Method),position=position_dodge(0.8))+theme_light()+
    scale_x_discrete(name ="Scale",limits=c("Leaf","Canopy","Plot"))+
    theme(axis.text.x = element_text(color="black",size=18),
          axis.text.y = element_text(color="black",size=18),axis.title=element_text(size=22,face="bold"),
          legend.text=element_text(size=16), legend.title=element_text(size=14,face = "bold"),
          legend.position=c(0.1,0.85))+
    scale_fill_grey(start = 0, end = 0.99)+labs(y="RMSEP")+ ylim(0.01,0.75)+
    annotate("text", x=2, y=0.74, label= "(c)",size=10)
  
  
}