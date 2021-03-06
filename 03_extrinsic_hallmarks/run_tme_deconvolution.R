library(nnet)
### Select train and test sets, dividing into 2/3 for train, 1/3 for test:
set.seed(19875)

folder_analysis<-getwd()

setwd('..')
current_dir<-getwd()
input_dir<-paste(current_dir,"/data",sep="")
output_dir<-paste(current_dir,"/output_dir",sep="")

setwd(input_dir)
load("PancancerTME.ssgsea.InputMultinom.RData")

input_mnlogit<-matrix_emt_and_tme[,c(1,3,4,6:ncol(matrix_emt_and_tme))]

randTrainingSet <- floor(0.70 * nrow(input_mnlogit))
train_ind <- sample(seq_len(nrow(input_mnlogit)), size = randTrainingSet)

dat.train <- input_mnlogit[train_ind,]
dat.test <- input_mnlogit[-train_ind,]

# dat.train.x <- as.matrix(dat.train[,1:(ncol(dat.train))])
# dat.train.y <- dat.train$states
# 
# dat.test.x <- as.matrix(dat.test[,1:(ncol(dat.test))])
# dat.test.y <- dat.test$states

# Fit multinomial logistic regression model:
glm.fit=multinom(states~0+., data=dat.train[,-c(1:2)])
summary(glm.fit)
z <- summary(glm.fit)$coefficients/summary(glm.fit)$standard.error
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

pmultinom<-t(data.frame(p))
RRmultinom<-t(summary(glm.fit)$coefficients)
colnames(RRmultinom)<-paste("RR",colnames(RRmultinom),sep="_")
out_multinom<-cbind(pmultinom,RRmultinom)

setwd(output_dir)
write.table(out_multinom,file="TME_multinomial_res_hmm_3states.ssgsea.txt",sep="\t",row.names=T,quote=F)


## Predictions:
predicted=predict(glm.fit,dat.test[,-c(1:2)],type="probs")
bpp=cbind(dat.test[,-c(colnames(dat.test)%in%c("samples"))], predicted)

bpp2 <- melt(bpp,id.vars = setdiff(colnames(bpp),c("epi","pEMT","mes")))

colnames(bpp2)[22:23] <- c("Prediction","Probability")


all_tme<-setdiff(colnames(bpp2)[-1],c("states","Prediction","Probability"))

library("gridExtra")

input_for_heatmap_tme<-data.frame()

pdf("TME_multinom_with_emt_ssgsea.pdf",width=10)

for(tmemn in 1:length(all_tme)){

current_tme<-all_tme[tmemn]
print(current_tme)

bpp3<-bpp2[,colnames(bpp2)%in%c(current_tme,"score_emt","Probability","Prediction")]
colnames(bpp3)[2]<-"tme"

input_for_heatmap_tme<-rbind(input_for_heatmap_tme,data.frame(immcell=current_tme,bpp3))

ptme <- ggplot(bpp3, aes(x = tme, y = Probability, colour = Prediction))+geom_point(shape=18,size=2)+scale_color_manual(values=c('#56b4e9','#E69f00','#EE0C0C')) + facet_grid(Prediction ~ ., scales="free")+ylim(c(0,1))+ylab("Probability")+labs(title=all_tme[tmemn]) + geom_smooth(method = "lm", formula = y ~ poly(x, 10),se=TRUE, linetype="dashed",color="blue")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

        
print(ptme)

}

dev.off()

pemt_tme<-rownames(out_multinom[out_multinom[,"pEMT"]<=0.01 & out_multinom[,"RR_pEMT"]>log(100),])
mes_tme<-rownames(out_multinom[out_multinom[,"mes"]<=0.01 & out_multinom[,"RR_mes"]>log(100),])

terms_to_consider<-unique(c(mes_tme,pemt_tme))

library(VennDetail)

ven <- venndetail(list(pemt = pemt_tme, mes = mes_tme))

write.table(file="TME_multinomial_res_hmm_3states.ssgsea.TERMSFORCHARTS.txt",result(ven, wide = TRUE),row.names=T,quote=F,sep="\t")

pdf("heatmap_multnom_tme.pdf",height=4)
p <- ggplot(input_for_heatmap_tme[input_for_heatmap_tme$immcell%in%terms_to_consider,],aes(x=tme,y=immcell,fill=Probability))+
  geom_tile()+ geom_vline(xintercept = 0, linetype="dotted",color = "black", size=0)+scale_fill_distiller(palette = "RdYlGn", direction = -1)+facet_grid(.~Prediction,scales="free_x")+theme_bw()
print(p)
dev.off()

pdf("heatmap_multnom_tme_scatterplotgrid.pdf",height=4)
p2<-ggplot(input_for_heatmap_tme[input_for_heatmap_tme$immcell%in%terms_to_consider,], aes(x=tme, y=Probability, color=Prediction)) +geom_point(shape=18,size=0.5)+ylim(c(0,1)) +geom_smooth(method = "lm", formula = y ~ poly(x, 10),se=TRUE, linetype="dashed",color="blue") + facet_grid(immcell~Prediction)
print(p2)
dev.off()

library(ggpubr)

input_for_boxplot_tme<-melt(input_mnlogit[,3:ncol(input_mnlogit)])
my_comparisons<-list(c("mes","epi"),c("epi","pEMT"),c("mes","pEMT"))

pdf("boxplot_TME_HMM_3_bytumor.ssgsea.selectedTME_multinom.pdf")

p1<-ggviolin(input_for_boxplot_tme[input_for_boxplot_tme$variable %in% terms_to_consider,], x = "states", y = "value", fill = "states",
         palette = c("#56b4e9", "#E69f00", "#EE0C0C"),
         add = "boxplot", add.params = list(fill = "white"),ylim=c(-0.7,1.1))+
         stat_compare_means(comparisons = my_comparisons,method="wilcox" ,label = "p.signif")+facet_wrap(~variable,nrow=5,ncol=5)+
         stat_summary(fun.y=median, geom="line", aes(group=1))

print(p1)

dev.off()

