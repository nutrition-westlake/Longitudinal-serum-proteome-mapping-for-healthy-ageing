library(glmnet)
library(randomForest)
library(vegan)
library(pROC)
library(sampling)
library(dplyr)
library(ggplot2)
library(splines)


########## LASSO model ##########
phase1 <- read.csv("phase1.csv",header=TRUE)
phase2 <- read.csv("phase2.csv",header=TRUE)
fh <- read.csv("fh.csv",header=TRUE)                                                                                   


set.seed(123)
phase1_train <- phase1 %>% group_by(id) %>% sample_n(size = 1)
phase1_test <- phase1[!phase1$idtime %in% phase1_train$idtime,]

xtrain=as.matrix(phase1_train[,4:90])
ytrain=as.matrix(phase1_train[,3])
xtestphase1 <- as.matrix(phase1_test[,4:90])
xtestphase1all<<- as.matrix(phase1[,4:90]) 
xtestphase2<- as.matrix(phase2[,4:90])
xtestfh<- as.matrix(fh[,4:90])

fit = glmnet(xtrain, ytrain, family="gaussian", nlambda=100, alpha=1)
print(fit)
plot(fit, xvar="lambda", label=FALSE)
fitcsv<-cbind.data.frame(fit$df,fit$dev.ratio,fit$lambda)
write.csv(fitcsv,file="fit.csv")

cvfit=cv.glmnet(xtrain,ytrain,nfolds=10)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
print(cvfit)
cvfit

ypredictphase1<-data.frame(phase1_test[,c(1,3)]) 
ypredictphase1all<-data.frame(phase1[,c(1,2,3)]) 
ypredictphase2<-data.frame(phase2[,c(1,3)]) 
ypredictfh<-data.frame(fh[,c(1,3)]) 

coefpro=as.data.frame(1:88,)
colnames(coefpro)[1]<-"logpro"

lambdalist<-read.csv("fit.csv")

for (i in 1:31) {
  x=lambdalist[i,4]
  ypredictphase1now<- predict(cvfit,newx=xtestphase1, s=x)
  ypredictphase2now<- predict(cvfit,newx=xtestphase2, s=x)
  ypredictfhnow<- predict(cvfit,newx=xtestfh, s=x)
  ypredictphase1<-cbind(ypredictphase1,ypredictphase1now) 
  ypredictphase2<-cbind(ypredictphase2,ypredictphase2now)
  ypredictfh<-cbind(ypredictfh,ypredictfhnow)
  
  ypredictcoefnow<-predict(cvfit, type="coefficients",newx=xtestphase1, s=x)
  coef<-as.matrix(ypredictcoefnow) 
  pro<-row.names(coef) 
  coefpro<-cbind(coefpro,coef)  
}

write.csv(coefpro,file="coefpro.csv")

###model performance
rypredictphase1<-data.frame(matrix(nrow=31,ncol=1))
colnames(rypredictphase1)[1]<-"rphase1"
rypredictphase2<-rypredictphase1
colnames(rypredictphase2)[1]<-"rphase2"
rypredictfh<-rypredictphase1
colnames(rypredictfh)[1]<-"rfh"

for (i in 3:33) {
  r=cor(ypredictphase1[,2],ypredictphase1[,i],method="pearson")
  rypredictphase1[i-2,1]<-r
}

for (i in 3:33) {
  r=cor(ypredictphase2[,2],ypredictphase2[,i],method="pearson")
  rypredictphase2[i-2,1]<-r
}

for (i in 3:33) {
  r=cor(ypredictfh[,2],ypredictfh[,i],method="pearson")
  rypredictfh[i-2,1]<-r
}

rypredict<-cbind(rypredictphase1,rypredictphase2,rypredictfh)
write.csv(rypredict,"rypredict.csv")


###calculate MAE
dypredictphase1<-data.frame(matrix(nrow=31,ncol=1))
colnames(dypredictphase1)[1]<-"mae"
dypredictphase2<-dypredictphase1
colnames(dypredictphase2)[1]<-"mae"
dypredictfh<-dypredictphase1
colnames(dypredictfh)[1]<-"mae"


for (i in 3:33) {
  ypredictphase1[,i+31]=ypredictphase1[,i]-ypredictphase1[,2]
  ypredictphase1[,i+31]=abs(ypredictphase1[,i+31])
  dypredictphase1[i-2,1]<-mean(ypredictphase1[,i+31])
}


for (i in 3:33) {
  ypredictphase2[,i+31]=ypredictphase2[,i]-ypredictphase2[,2]
  ypredictphase2[,i+31]=abs(ypredictphase2[,i+31])
  dypredictphase2[i-2,1]<-mean(ypredictphase2[,i+31])
}


for (i in 3:33) {
  ypredictfh[,i+31]=ypredictfh[,i]-ypredictfh[,2]
  ypredictfh[,i+31]=abs(ypredictfh[,i+31])
  dypredictfh[i-2,1]<-mean(ypredictfh[,i+31])
}

dypredict<-cbind(dypredictphase1,dypredictphase2,dypredictfh)
write.csv(dypredict,"dypredict.csv")

write.csv(ypredictphase1,"ypredictphase1.csv")
write.csv(ypredictphase2,"ypredictphase2.csv")
write.csv(ypredictfh,"ypredictfh.csv")
ypredictphase1allnow<- predict(cvfit,newx=xtestphase1all, s=lambdalist[31,4])
ypredictphase1all<-cbind(ypredictphase1all,ypredictphase1allnow)
write.csv(ypredictphase1all,"1ypredictphase1all.csv")



########## random forest model ##########

phase12 <- read.csv("phase12_for_ml.csv",header=TRUE)
FH <- read.csv("FH_for_ml.csv",header=TRUE)
phase12$healthy<-as.factor(phase12$healthy)
phase12$sex<-as.factor(phase12$sex)

phase1<-subset(phase12,phase==1)
phase1_train<-subset(phase1,time==1)
phase1_test<-subset(phase1,time>1)

phase2<-subset(phase12,phase==2)
phase2_test<-subset(phase2,time==1)
phase2_test23<-subset(phase2,time>1)


###random forest models using all proteins
train_allproteins <-data.frame(phase1_train[,c(6,10:417)])
test_allproteins <-data.frame(phase2_test[,c(6,10:417)])

set.seed(10)
rf_allproteins<-randomForest(healthy~.,data=train_allproteins,importance=TRUE,ntree=1000,na.action=na.omit)
rf_allproteins

predict_train_allproteins <- predict(rf_allproteins,newdata=train_allproteins,type="prob")
predictc_train_allproteins <- predict(rf_allproteins,newdata=train_allproteins)

predict_test_allproteins <- predict(rf_allproteins,newdata=test_allproteins,type="prob")
predictc_test_allproteins<- predict(rf_allproteins,newdata=test_allproteins)

table(train_allproteins$healthy,predictc_train_allproteins,dnn=c("Actual","Predicted"))
table(test_allproteins$healthy,predictc_test_allproteins,dnn=c("Actual","Predicted"))

roc_allproteins<-roc(test_allproteins$healthy,as.numeric(predict_test_allproteins[,2]))
plot.roc(test_allproteins$healthy,as.numeric(predict_test_allproteins[,2]),print.auc=TRUE,ci=F,col="#d33524")


###random forest models using ageing-related proteins
train_ageingproteins <-data.frame(phase1_train[,c(6,10:95)])
test_ageingproteins <-data.frame(phase2_test[,c(6,10:95)])

set.seed(10)
rf_ageingproteins<-randomForest(healthy~.,data=train_ageingproteins,mtry=4,importance=TRUE,ntree=1000,na.action=na.omit)
rf_ageingproteins
plot(rf_ageingproteins)
varImpPlot(rf_ageingproteins, main = NULL)

rf_importance<-data.frame(importance(rf_ageingproteins))
head(rf_importance)
rf_importance <- rf_importance[order(rf_importance$MeanDecreaseAccuracy, decreasing = TRUE), ]
rf_importance1<- rf_importance[order(rf_importance$MeanDecreaseGini, decreasing = TRUE), ]
write.csv(rf_importance,file="rf_importance.csv")

predict_train_ageingproteins <- predict(rf_ageingproteins,newdata=train_ageingproteins,type="prob")
predictc_train_ageingproteins <- predict(rf_ageingproteins,newdata=train_ageingproteins)
predict_test_ageingproteins <- predict(rf_ageingproteins,newdata=test_ageingproteins,type="prob")
predictc_test_ageingproteins<- predict(rf_ageingproteins,newdata=test_ageingproteins)

table(train_ageingproteins$healthy,predictc_train_ageingproteins,dnn=c("Actual","Predicted"))
table(test_ageingproteins$healthy,predictc_test_ageingproteins,dnn=c("Actual","Predicted"))

roc_ageingproteins<-roc(test_ageingproteins$healthy,as.numeric(predict_test_ageingproteins[,2]))
plot.roc(test_ageingproteins$healthy,as.numeric(predict_test_ageingproteins[,1]),print.auc=TRUE,ci=F,col="#d33524")


###cross-validation
rf_cv<-replicate(5,rfcv(train_ageingproteins[2:87],train_ageingproteins$healthy,cv.fold=10,step=1.1),simplify=FALSE)
rf_cv

rf_cv<-data.frame(sapply(rf_cv,'[[', 'error.cv'))
rf_cv$pros<-rownames(rf_cv)
rf_cv<-reshape2::melt(rf_cv,id='pros')
rf_cv$pros<-as.numeric(as.character(rf_cv$pros))

p<- ggplot(rf_cv, aes(pros, value)) +
  geom_smooth(se = TRUE,method = 'glm', color="#fd9200",formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of variables', y = 'Cross-validation error')
p
p + geom_vline(xintercept=22,color="darkgrey",linetype="dashed",size=0.5)

selectvar<-row.names(rf_importance)[1:22]
train_ageingproteins_top<-train_ageingproteins[,c(selectvar,'healthy')]
test_ageingproteins_top<-test_ageingproteins[,c(selectvar,'healthy')]

set.seed(10)
rf_top<-randomForest(healthy~.,data=train_ageingproteins_top,importance=TRUE,mtry=3,ntree=1000,na.action=na.omit)

predict_train_ageingproteins_top <- predict(rf_top,newdata=train_ageingproteins_top,type="prob")
predictc_train_ageingproteins_top <- predict(rf_top,newdata=train_ageingproteins_top)
predict_test_ageingproteins_top <- predict(rf_top,newdata=test_ageingproteins_top,type="prob")
predictc_test_ageingproteins_top <- predict(rf_top,newdata=test_ageingproteins_top)

table(train_ageingproteins_top$healthy,predictc_train_ageingproteins_top,dnn=c("Actual","Predicted"))
table(test_ageingproteins_top$healthy,predictc_test_ageingproteins_top,dnn=c("Actual","Predicted"))

roc_ageingproteins_top<-roc(test_ageingproteins_top$healthy,as.numeric(predict_test_ageingproteins_top[,2]))
plot.roc(test_ageingproteins_top$healthy,as.numeric(predict_test_ageingproteins_top[,2]),print.auc=TRUE,ci=F,col="#d33524")

plot(roc_allproteins,col="#fd9200",grid=c(0.2, 0.2))
plot(roc_ageingproteins,add=TRUE,col="#d33524")
plot(roc_ageingproteins_top,add=TRUE,col="#1576ae")
legend("bottomright", legend=c("AUC=0.72(409 proteins)","AUC=0.70 (86 ageing proteins)","AUC=0.70 (22 ageing proteins)"),col=c("#fd9200","#d33524","#1576ae"),lty=1)

###generate PHAS
phase1_top<-phase1[,c('healthy',selectvar)]
phase2_top<-phase2[,c('healthy',selectvar)]

predict_phase1 <- predict(rf_top,newdata=phase1_top,type="prob")
predict_phase1c <- predict(rf_top,newdata=phase1_top)

predict_phase2 <- predict(rf_top,newdata=phase2_top,type="prob")
predict_phase2_c <- predict(rf_top,newdata=phase2_top)


ypredict_phase1 = data.frame(prob=predict_phase1[,2],phase1_top[,1],phase1)
write.csv(ypredict_phase1,file="phase1_prob.csv")

ypredict_phase2 = data.frame(prob=predict_phase2[,2],phase2_top[,1],phase2)
write.csv(ypredict_phase2,file="phase2_prob.csv")


###random forest models using age sex bmi and 22 proteins
train_full<-data.frame(phase1_train[,c("healthy","ageadj","bmiadj","sexadj",selectvar)])
test_full<-data.frame(phase2_test[,c("healthy","ageadj","bmiadj","sexadj",selectvar)])

set.seed(10)
rf_full<-randomForest(healthy~.,data=train_full,importance=TRUE,ntree=1000,na.action=na.omit)

predict_train_full <- predict(rf_full,newdata=train_full,type="prob")
predictc_train_full <- predict(rf_full,newdata=train_full)
predict_test_full <- predict(rf_full,newdata=test_full,type="prob")
predictc_test_full <- predict(rf_full,newdata=test_full)

table(train_full$healthy,predictc_train_full,dnn=c("Actual","Predicted"))
table(test_full$healthy,predictc_test_full,dnn=c("Actual","Predicted"))

roc_full<-roc(test_full$healthy,as.numeric(predict_test_full[,2]))
plot.roc(test_full$healthy,as.numeric(predict_test_full[,2]),print.auc=TRUE,ci=F,col="#d33524")


###random forest models using age sex bmi
train_agesexbmi <-data.frame(phase1_train[,c(6:9)])
test_agesexbmi <-data.frame(phase2_test[,c(6:9)])

set.seed(10)
rf_agesexbmi<-randomForest(healthy~.,data=train_agesexbmi,importance=TRUE,ntree=1000,na.action=na.omit)
rf_agesexbmi

predict_train_agesexbmi <- predict(rf_agesexbmi,newdata=train_agesexbmi,type="prob")
predictc_train_agesexbmi <- predict(rf_agesexbmi,newdata=train_agesexbmi)
predict_test_agesexbmi <- predict(rf_agesexbmi,newdata=test_agesexbmi,type="prob")
predictc_test_agesexbmi<- predict(rf_agesexbmi,newdata=test_agesexbmi)

table(train_agesexbmi$healthy,predictc_train_agesexbmi,dnn=c("Actual","Predicted"))
table(test_agesexbmi$healthy,predictc_test_agesexbmi,dnn=c("Actual","Predicted"))

roc_agesexbmi<-roc(test_agesexbmi$healthy,as.numeric(predict_test_agesexbmi[,2]))
plot.roc(test_agesexbmi$healthy,as.numeric(predict_test_agesexbmi[,2]),print.auc=TRUE,ci=F,col="#d33524")

plot(roc_full,col="#1576ae",grid=c(0.2, 0.2))
plot(roc_ageingproteins_top,add=TRUE,col="#d33524")
plot(roc_agesexbmi,add=TRUE,col="#fd9200")
legend("bottomright", legend=c("AUC=0.72 (Age,sex,BMI+22 ageing proteins)","AUC=0.70 (PHAS/22 ageing proteins)","AUC=0.62 (Age,sex,BMI)"),col=c("#1576ae","#d33524","#fd9200"),lty=1)


predict_FH<-predict(rf_top,newdata=FH,type="prob")
ypredict_FH = data.frame(prob=predict_FH[,2],FH)
write.csv(ypredict_FH,file="FH_prob.csv")



########## PERMANOVA for 22 proteins##########
prodist<-read.csv("prodist.csv",head=T,row.names=1)
logprodist<-read.csv("logprodist.csv",head=T,row.names=1)
stdlogprodist<-read.csv("stdlogprodist.csv",head=T,row.names=1)
life<-read.csv("life.csv",head=T,row.names=1)
life$smokeadj<-as.factor(life$smokeadj)
life$sexadj<-as.factor(life$sexadj)
life$alcadj<-as.factor(life$alcadj)
diet<-read.csv("diet.csv",header = TRUE,row.names = 1)
metagenomic<-read.csv("metagenomic_219species.csv",header=TRUE,row.names = 1)
genome<-read.csv("genome.csv",header=TRUE,row.names = 1)

###intrinsic and lifestyle
set.seed(10)
div_demo<-adonis2(prodist ~ ageadj+sexadj+bmiadj, data = life,distance = "bray",permutations = 999,by=NULL)
div_demo<-adonis2(prodist ~ ageadj+sexadj+bmiadj, data = life,distance = "bray",permutations = 999)
div_demo

set.seed(10)
div_life<-adonis2(prodist ~ smokeadj+metadj+alcadj, data = life,distance = "bray",permutations = 999,by=NULL)
div_life<-adonis2(prodist ~ smokeadj+metadj, data = life,distance = "bray",permutations = 999)
div_life<-adonis2(prodist ~ smokeadj+metadj, data = life,distance = "bray",permutations = 999,by=NULL)
div_life


###diet
result_diet=as.data.frame(matrix(nrow=15,ncol=4)) ####璁剧疆寰幆缁撴灉p鍊兼暟鎹锛?219涓洪暱搴?
colnames(result_diet)[1]<-"SumOfSqs"
colnames(result_diet)[2]<-"R2"
colnames(result_diet)[3]<-"F"
colnames(result_diet)[4]<-"pvalue"

for (i in 1:15) {
  x=colnames(diet)[i]
  form<-as.formula(paste("prodist",x,sep="~"))
  set.seed(10)
  div_diet_individual<-adonis2(form, data = diet,distance = "bray",permutations = 999)
  result_diet[i,1]<-div_diet_individual$SumOfSqs[1]
  result_diet[i,2]<-div_diet_individual$R2[1]
  result_diet[i,3]<-div_diet_individual$F[1]
  result_diet[i,4]<-div_diet_individual$Pr[1]
}

div_diet_step2<-adonis2(prodist ~ vegetables_withoutpotato_g+juice_g+fish_g, data = diet,distance = "bray",permutations = 999)
div_diet_step2
div_diet_step3<-adonis2(prodist ~ vegetables_withoutpotato_g+juice_g, data = diet,distance = "bray",permutations = 999)
div_diet_step3
div_diet_final<-adonis2(prodist ~ vegetables_withoutpotato_g+juice_g, data = diet,distance = "bray",permutations = 999,by=NULL)
div_diet_final


###microbiota
result_species=as.data.frame(matrix(nrow=219,ncol=4)) ####璁剧疆寰幆缁撴灉p鍊兼暟鎹锛?219涓洪暱搴?
colnames(result_species)[1]<-"SumOfSqs"
colnames(result_species)[2]<-"R2"
colnames(result_species)[3]<-"F"
colnames(result_species)[4]<-"pvalue"

for (i in 1:219) {
  x=colnames(metagenomic)[i]
  form<-as.formula(paste("prodist",x,sep="~"))
  set.seed(10)
  div_metagenomics<-adonis2(form, data = metagenomic,distance = "bray",permutations = 999)
  result_species[i,1]<-div_metagenomics$SumOfSqs[1]
  result_species[i,2]<-div_metagenomics$R2[1]
  result_species[i,3]<-div_metagenomics$F[1]
  result_species[i,4]<-div_metagenomics$Pr[1]
}

result_species$qvalue<-p.adjust(result_species$pvalue,method="BH")
result_species_step2<-subset(result_species,pvalue < 0.05)

#step 2
set.seed(10)
div_species_step2<-adonis2(prodist~v2+v24+v32+v35+v43+v59+v92+v100+v103+v132+v143+v147+v153+v154+v155+v161+v165+v167+v171+v183+v209+v219, data = metagenomic,distance = "bray",permutations = 999)
div_species_step2

#step 3
set.seed(10)
div_species_step3<-adonis2(prodist~v2+v24+v32+v35+v43+v59+v92+v100+v103+v143+v147+v155+v161+v209+v219, data = metagenomic,distance = "bray",permutations = 999)
div_species_step3

#step 4
set.seed(10)
div_species_final<-adonis2(prodist~v2+v24+v32+v35+v43+v59+v92+v100+v103+v143+v147+v155+v161+v209+v219, data = metagenomic,distance = "bray",permutations = 999,by=NULL)
div_species_final


###genetic
result_genome=as.data.frame(matrix(nrow=65,ncol=4))
colnames(result_genome)[1]<-"SumOfSqs"
colnames(result_genome)[2]<-"R2"
colnames(result_genome)[3]<-"F"
colnames(result_genome)[4]<-"pvalue"

for (i in 1:65) {
  x=colnames(genome)[i]
  form<-as.formula(paste("prodist",x,sep="~"))
  set.seed(10)
  div_genome<-adonis2(form, data = genome,distance = "bray",permutations = 999)
  result_genome[i,1]<-div_genome$SumOfSqs[1]
  result_genome[i,2]<-div_genome$R2[1]
  result_genome[i,3]<-div_genome$F[1]
  result_genome[i,4]<-div_genome$Pr[1]
}

result_genome_step2<-subset(result_genome,pvalue < 0.05)

####step 2
set.seed(10)
div_genome_step2<-adonis2(prodist~snp1+snp2+snp7+snp8+snp9+snp10+snp14+snp15+snp20+snp21+snp22+snp23+snp24+snp25+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp38+snp43+snp55+snp56+snp59+snp60+snp61+snp62+snp63+snp64+snp65, data = genome,distance = "bray",permutations = 999)
div_genome_step2

####step 3
set.seed(10)
div_genome_step3<-adonis2(prodist~snp1+snp7+snp14+snp20+snp21+snp22+snp23+snp25+snp27+snp28+snp43+snp59+snp61+sn3, data = genome,distance = "bray",permutations = 999)
div_genome_step3

####step 4
set.seed(10)
div_genome_step4<-adonis2(prodist~snp1+snp7+snp14+snp20+snp21+snp22+snp23+snp25+snp27+snp28+snp43+snp59+snp61+snp63, data = genome,distance = "bray",permutations = 999,by=NULL)
div_genome_step4

####step 4 
set.seed(10)
div_genome_step4_log<-adonis2(logprodist~snp1+snp7+snp14+snp20+snp21+snp22+snp23+snp25+snp27+snp28+snp43+snp59+snp61+snp63, data = genome,distance = "bray",permutations = 999)
div_genome_step4_log


########## Variance for PHAS ##########
PHAS <- read.csv("PHAS.csv",header=TRUE,row.names = 1)
y<-as.matrix(PHAS[,1])
x_life<-as.matrix(life)
x_diet<-as.matrix(diet)
x_metagenomic<-as.matrix(metagenomic)
x_genome<-as.matrix(genome)
data<-cbind(PHAS,ypredict,life,diet,metagenomic,genome)
write.csv(data,"multiomics.csv",row.names = TRUE)

PHAS_life<-cbind(PHAS,life)
#intrinsic
set.seed(123)
cvfit=cv.glmnet(demo,y,nfolds=10)
cvfit
coef_demo_HAS<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact= F)
coef_demo_HAS
lm<-lm(formula = prob ~  ageadj+sexadj+bmiadj, data = PHAS_life)
summary(lm)

#lifestyle
set.seed(123)
cvfit=cv.glmnet(lifestyle,y,nfolds=10)
coef_life_HAS<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact= F)
coef_life_HAS
cvfit
lm<-lm(formula = prob ~  metadj, data = PHAS_life)
summary(lm)

#diet
set.seed(123)
cvfit=cv.glmnet(x_diet,y,nfolds=10)
print(cvfit)
coef_diet<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact= F)
coef_diet1<-as.matrix(coef_diet)
coef_diet1

PHAS_diet<-cbind(PHAS,x_diet)
lm<-lm(formula = prob ~  vegetables_withoutpotato_g+coffee+teaadj, data = PHAS_diet)
summary(lm)

#metagenomic
set.seed(123)
cvfit=cv.glmnet(x_metagenomic,y,nfolds=10)
cvfit
coef_metagenomic<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact= F)
coef_metagenomic<-as.matrix(coef_metagenomic)
coef_metagenomic<-as.data.frame(coef_metagenomic)
coef_metagenomic=cbind(var=row.names(coef_metagenomic), coef_metagenomic)
coef_metagenomic<-coef_metagenomic[-1,]
colnames(coef_metagenomic)[2]<-"coef"
coef_metagenomic_subset<-subset(coef_metagenomic,coef!=0)
coef_metagenomic_subset
var<-paste(rownames(coef_metagenomic_subset),collapse = "+")
form<-as.formula(paste(colnames(PHAS_metagenomic)[1],var,sep="~"))
lm<-lm(formula = form, data = PHAS_metagenomic)
summary(lm)
form

#genetic
set.seed(123)
cvfit=cv.glmnet(x_genome,y,nfolds=10)
cvfit
coef_genome<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact= F)
coef_genome<-as.matrix(coef_genome)
coef_genome<-as.data.frame(coef_genome)
coef_genome=cbind(var=row.names(coef_genome), coef_genome)
coef_genome<-coef_genome[-1,]
colnames(coef_genome)[2]<-"coef"
coef_genome_subset<-subset(coef_genome,coef!=0)
coef_genome_subset
var<-paste(rownames(coef_genome_subset),collapse = "+")
form<-as.formula(paste(colnames(PHAS_genome)[1],var,sep="~"))
lm<-lm(formula = form, data = PHAS_genome)
summary(lm)
form
