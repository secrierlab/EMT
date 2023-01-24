petalChartGMT<-function(X=CCLE_pseudotime_EMT_metpot_HMM,var_group="hmm_states",title_string,output){
  
  RadarTheme<-theme(panel.background=element_blank(),
                    legend.position="bottom",legend.title=element_blank(),legend.direction="vertical",
                    axis.text.x = element_blank(),
                    axis.line.x=element_line(size=0.5),
                    panel.grid.major=element_line(size=0.3,linetype = 2,colour="grey"))
  
  X2<-X[,which(colnames(X)%in%c(var_group,"penetrance"))]
  X2_penetrance<-aggregate(.~hmm_states,X2,mean)
  width_petal<-X2_penetrance[,2]
  
  X2<-X[,which(colnames(X)%in%c(var_group,"mean"))]
  X2_metpot<-aggregate(.~hmm_states,X2,mean)
  mean_metpot<-X2_mean[,2]
  
  p1<-ggplot(X, aes(x = hmm_states, y = mean, fill = hmm_states))+scale_fill_manual(values = alpha(c("red2", "steelblue2","orange2"), 0.5))+ geom_hline(yintercept=-2, linetype="dashed", color = "red")+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    geom_hline(yintercept=-2, linetype="dashed", color = "red")+
    stat_summary(geom = "bar", fun = "max", position = position_dodge(0.9),width= width_petal)+
    stat_summary(geom = "bar", fun = "min", position = position_dodge(0.9),width= width_petal)+
    stat_summary(
      mapping = aes(x = hmm_states, y = mean),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, color="black",size = 0.1)+
    coord_polar()+RadarTheme+ggtitle(title_string)
  
  p2<-ggplot(X, aes(x = hmm_states, y = mean, fill = hmm_states), colour="black")+scale_fill_manual(values = alpha(c("red2", "steelblue2","orange2"), 0.5))+scale_colour_manual("black")+ geom_hline(yintercept=-2, linetype="dashed", color = "red")+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    geom_hline(yintercept=-2, linetype="dashed", color = "red")+
    stat_summary(geom = "bar", fun = "max", position = position_dodge(0.9),width= width_petal)+
    stat_summary(geom = "bar", fun = "min", position = position_dodge(0.9),width= width_petal)+
    stat_summary(
      mapping = aes(x = hmm_states, y = mean),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, color="black",size = 0.1)+facet_wrap(~target_tissue,nrow=3,ncol=2)+coord_polar()+RadarTheme+ggtitle(title_string)
  
  pdf(output)
  
  print(p1)
  
  print(p2)
  
  dev.off()
  
}


petalChartGMT_EMT_for_METPOT<-function(X=CCLE_pseudotime_EMT_metpot_HMM,var_group="status_metpot",title_string,output){
  
  RadarTheme<-theme(panel.background=element_blank(),
                    legend.position="bottom",legend.title=element_blank(),legend.direction="vertical",
                    axis.text.x = element_blank(),
                    axis.line.x=element_line(size=0.5),
                    panel.grid.major=element_line(size=0.3,linetype = 2,colour="grey"))
  
  X2<-X[,which(colnames(X)%in%c(var_group,"penetrance"))]
  X2_penetrance<-aggregate(.~status_metpot,X2,mean)
  width_petal<-X2_penetrance[,2]
  
  p1<-ggplot(X, aes(x = status_metpot, y = emt_score, fill = status_metpot))+scale_fill_manual(values = alpha(c( "steelblue2","orange2","red2"), 0.5))+ geom_hline(yintercept=-2, linetype="dashed", color = "red")+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    stat_summary(geom = "bar", fun = "max", position = position_dodge(0.9),width= width_petal)+
    stat_summary(geom = "bar", fun = "min", position = position_dodge(0.9),width= width_petal)+
    stat_summary(
      mapping = aes(x = status_metpot, y = emt_score),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, color="black",size = 0.1)+
    coord_polar()+RadarTheme+ggtitle(title_string)
  
  p2<-ggplot(X, aes(x = status_metpot, y = emt_score, fill = status_metpot), colour="black")+scale_fill_manual(values = alpha(c("steelblue2","orange2","red2"), 0.5))+scale_colour_manual("black")+ geom_hline(yintercept=-2, linetype="dashed", color = "red")+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    stat_summary(geom = "bar", fun = "max", position = position_dodge(0.9),width= width_petal)+
    stat_summary(geom = "bar", fun = "min", position = position_dodge(0.9),width= width_petal)+
    stat_summary(
      mapping = aes(x = status_metpot, y = emt_score),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, color="black",size = 0.1)+facet_wrap(~target_tissue,nrow=3,ncol=2)+coord_polar()+RadarTheme+ggtitle(title_string)
  
  pdf(output)
  
  print(p1)
  
  print(p2)
  
  dev.off()
  
}
