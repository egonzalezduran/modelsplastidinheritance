####################### R Code for statistical models (GLM) ########################################
## "Control of plastid inheritance by environmental and genetic factors"                         ###
## Chung, K.P ; Gonzalez-Duran, E; Ruf, S.; Endries, P. and Bock, R (2022)                       ###
## Publishing in Nature Plants                                                                   ###
## Version 17.11.2022 by Enrique Gonzalez-Duran,                                                 ###
## Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                     ###  
## R version 3.5.3 accessed through R studio  22.07.2022 Build 576                               ###   
################################################################################################ ###   

#### List of statistical models as listed in Extended Data Tables 1 and 2 and referred by in the manuscript


## Model 1: Effect of stres treatments during pollen development over paternal plastid transmission
## Model 2: Effect of chilling stress and stage of visualization over inclusion of plastids in GCs 
## Model 3: Effect of chilling stress and dpd1 genotype over paternal plastid transmission across experiments
## Model 4: Effect of chilling stress on pollen viability
## Model 5: Effect of the dpd1 mutation on pollen viability


#### this code, run from top to bottom, should lead to the 5 plots that were used to assemble figures, and all calculations listed in the manuscript
##This code is set to work in R version 3.5.3. It was tested in R v.4.1.2 with minimal changes in the plotting behavior

##List of necessary additional packages:

# MASS        v. 7.3-51.1
# dplyr       v. 0.8.0.1
# ggplot2     v. 3.1.1
# ciTools     v. 0.5.1
# insight     v. 0.2.0
# pscl        v. 1.5.2
# MuMIN       v. 1.43.6
# rstudioapi  v. 0.10

wants <- c("MASS","dplyr", "ggplot2","insight","pscl","ciTools", "ciTools","MuMIn","rstudioapi") ## searches for and installs the packages
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

library(MASS) 
library(dplyr) 
library(ggplot2)
library(ciTools)
library(pscl)
library(insight)
library(MuMIn)
library(rstudioapi)

## List of necessary data files:
## data_for_model_1.txt
## data_for_model_2.txt
## data_for_model_3.txt
## data_for_model_4.txt
## data_for_model_5.txt


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets the location of the file as the working directory, works in R Studio
getwd()


########################################################################## Model 1 ############################################################################
## Necessary data file: data_for_model_1.txt

#Datapoint Treatment Total.seedlings Resistant 
#A Greenhouse 942992 16 
#B Greenhouse 238423 6
#C Greenhouse 637379 10
#D Greenhouse 275794 1
#E High.light 152719 2
#F High.light 146631 1
#G High.light 174965 0
#H Heat 247451 0
#I Heat 248615 1
#J Heat 231540 3
#K Drought 173077 2
#L Drought 131329 4
#M Drought 167625 4
#N Cold 167169 286
#O Cold 147661 402
#P Cold 164400 463

large.stress.exp <- read.table("data_for_model_1.txt", header = TRUE) #passes table from txt file to a data frame
LSexp1data <- as.data.frame(large.stress.exp)
LSexp1df <- as.data.frame(large.stress.exp)
LSexp1data$Treatment = factor(LSexp1df$Treatment, levels=unique(LSexp1df$Treatment))

### in the following glm models, parameter effects are provided as natural logarithms of the effect, they are converted to base 10 later on


LSexp1.nb.M1 <- glm.nb(Resistant ~ Treatment + offset(log(Total.seedlings)),
                contrasts=list(Treatment=contr.treatment(n=5,base=1)), data=LSexp1data)


LSexp1.poisson.m1.2 <- glm(Resistant ~ Treatment  + offset(log(Total.seedlings)), family= poisson(link="log"), 
               contrasts=list(Treatment=contr.treatment(5,base=1)), data=LSexp1data ) 

#Negative binomial model plotting
##M1 parameters

LSexp1df$frequency.per.datapoint <- LSexp1df$Resistant / LSexp1df$Total.seedlings ### this expresses frequencies per data point, in natural log
LSexp1df$l10frequency.per.datapoint <- sapply(LSexp1df$freq, function (x) {log(x,base = 10)}) ### this expresses frequencies per datapoint in log 10
gpM1<-get_parameters(LSexp1.nb.M1) ## this function extracts the parameter values from the model
paramver.M1 <-as.vector(gpM1[,2]) ## vector with the parameter values

##M1 confidence intervals (of parameters)

summary.glm(LSexp1.nb.M1)
confintM1<- as.data.frame(confint.default(LSexp1.nb.M1)) ## puts the confidence intervals of effect values to a dataframe. The effects (and therefore, the intervals) are relative to the intercept (greenhouse treatment)
confintM1adj<- data.frame(matrix(data=NA,nrow=5,ncol=2)) 
confintM1adj[1,] <- confintM1[1,]
confintM1adj[2,] <- (confintM1[1,] + confintM1[2,])
confintM1adj[3,] <- (confintM1[1,] + confintM1[3,])
confintM1adj[4,] <- (confintM1[1,] + confintM1[4,])
confintM1adj[5,] <- (confintM1[1,] + confintM1[5,]) ## this is the transformation of the CIs of the effects (relative to greenhouse) to CIs of the frequencies (fixed), expressed still in natural log
confintM1lower <-as.vector(confintM1adj[,1]) ## this puts both CIs ordered their own dataframe
confintM1upper <-as.vector(confintM1adj[,2]) ##


#Assembly of a table with the means fitted by the model (they are equal to the frequency taking all replicates of a condition together)
df.fitted.meansM1 <- as.data.frame(matrix(data=NA, nrow=1, ncol=16))
df.fitted.meansM1[1,] <- fitted(LSexp1.nb.M1)
df.fitted.meansM1<- rbind(df.fitted.meansM1, c(1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
df.fitted.meansM1<- rbind(df.fitted.meansM1, c(LSexp1df$Total.seedlings))
df.fitted.meansM1<- rbind(df.fitted.meansM1,(as.numeric(fitted(LSexp1.nb.M1)/LSexp1df$Total.seedlings)))
rownames(df.fitted.meansM1)<-c("fitted.mean","grp","Total.Session","fitted.frequency")
colnames(df.fitted.meansM1)<-c("G1","G2","gPM3","G4","HL1","HL2","HL3","HT1","HT2","HT3","D1","D2","paramver.M2","C1","C2","C3")
df.fitted.meansM1<-as.data.frame(t(df.fitted.meansM1))
df.fitted.meansM1 ### this table is used in the next dataframe, ultimately will help with plotting the graph
###

###dataframe with efect details. 
df.effects.M1 <- as.data.frame(matrix(data=NA, nrow=5, ncol=8), rownames = c("G","HL","HT","D","C"))
colnames(df.effects.M1) <- c("Treatment","lower","upper","mean.per.treatment", "l10mean.per.treatment","grp","fitted.frequency","fitted.from.param")
df.effects.M1[,1]  <-  c("Greenhouse","High.light","Heat","Drought","Cold")
df.effects.M1[,2]  <-  log(exp(confintM1lower), base=10)
df.effects.M1[,3]  <-  log(exp(confintM1upper), base=10)
df.effects.M1[,4] <- c(33/2094588,3/474315,4/727606,10/472031,1161/479230) 
df.effects.M1[,5] <- log(as.vector(df.effects.M1[,4]), base = 10)
df.effects.M1[,6] <- c(1,2,3,4,5)
df.effects.M1[1,7] <- df.fitted.meansM1[[1,"fitted.frequency"]] #### frequency of paternal transmission per treatment fitted by the model taken from the points. 
df.effects.M1[2,7] <- df.fitted.meansM1[[5,"fitted.frequency"]]
df.effects.M1[3,7] <- df.fitted.meansM1[[8,"fitted.frequency"]]
df.effects.M1[4,7] <- df.fitted.meansM1[[11,"fitted.frequency"]]
df.effects.M1[5,7] <- df.fitted.meansM1[[14,"fitted.frequency"]]
df.effects.M1[1,8] <- exp(paramver.M1[1]) ### each of the parameter values is transformed here from relative to the greenhouse to the frequency according to the model. this is an alternative way to calculate the means per treatment, which must be equal to the ones fitted by the model (column 7)
df.effects.M1[2,8] <- exp(paramver.M1[1] + paramver.M1[2]) 
df.effects.M1[3,8] <- exp(paramver.M1[1] + paramver.M1[3])
df.effects.M1[4,8] <- exp(paramver.M1[1] + paramver.M1[4])
df.effects.M1[5,8] <- exp(paramver.M1[1] + paramver.M1[5])


df.for.plotting.M1 <- as.data.frame(dplyr::bind_rows(LSexp1df, df.effects.M1)) 
df.for.plotting.M1$l10fitted.frequency <- as.numeric(sapply(df.for.plotting.M1$fitted.frequency, function (x) {log(x,base = 10)})) ## this is the mean assigned by the model expressed in log10, goes in the black horizontal bars
 ### single dataframe with all needed for the plot

####### plotting 
plotM1 <- ggplot(df.for.plotting.M1, aes(x=factor(Treatment,level=c("Greenhouse","High.light","Heat","Drought","Cold")), y=l10frequency.per.datapoint, color= Treatment)) +
  geom_point(aes(size = Total.seedlings), position = position_dodge2(width=0.8), alpha=0.6, show.legend= TRUE) + #plots the points according to the previous line
  scale_size_area(max_size = 8) +
  xlab(" ") +
  ylab("Rate of paternal plastid transmission (log10)") +
  geom_errorbar(aes(x= factor(Treatment, level=c("Greenhouse","High.light","Heat","Drought","Cold")), ymin = l10fitted.frequency, ymax = l10fitted.frequency), col="#000000", linetype=1, size=0.5) + #plots the means
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_ribbon(aes( ymin=lower, ymax=upper), linetype=1, show.legend= FALSE, size=1) +  #plots the CI
  theme(axis.line = element_line(colour = "black"))  +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(limits=c(-8.5,0), expand = c(0,0)) 
plotM1

########################################################################## Model 2 ############################################################################
#Necessary data file: data_for_model_2.txt

#Group	Try	Number.sessions	Total.Session	Experiment	Temperature	Total.pollen	Included	Not.Included
#1	K	4	18	EPM	20C	56	2	16
#1	L	4	6	EPM	20C	56	1	5
#1	M	4	7	EPM	20C	56	2	5
#1	N	4	25	EPM	20C	56	2	23
#2	O	4	14	EPM	10C	65	4	10
#2	P	4	26	EPM	10C	65	8	18
#2	Q	4	12	EPM	10C	65	4	8
#2	R	4	13	EPM	10C	65	4	9
#3	A	3	25	GPT	20C	61	3	22
#3	B	3	26	GPT	20C	61	7	19
#3	C	3	10	GPT	20C	61	1	9
#4	D	7	17	GPT	10C	119	5	12
#4	E	7	16	GPT	10C	119	5	11
#4	F	7	18	GPT	10C	119	8	10
#4	G	7	20	GPT	10C	119	6	14
#4	H	7	21	GPT	10C	119	13	8
#4	I	7	20	GPT	10C	119	14	6
#4	J	7	7	GPT	10C	119	5	2

plastid.in.gc <- read.table("data_for_model_2.txt", header = TRUE)
pGCexpdf <- as.data.frame(plastid.in.gc)
pGCexpdf$Rate <- pGCexpdf$Included / pGCexpdf$Total.Session



pGCbinomial.M2 <- glm(Rate ~ Temperature + Experiment, weights = Total.Session, family= binomial(link= "log"),
                      contrasts=list(Experiment=contr.treatment(2,base=1), Temperature=contr.treatment(2,base=2)), data=pGCexpdf)



pGCbinomial.m2.2 <- glm(Rate ~ Temperature + Experiment + Temperature:Experiment, weights = Total.Session, family= binomial(link= "log"),
                     contrasts=list(Experiment=contr.treatment(2,base=1), Temperature=contr.treatment), data=pGCexpdf)



Inc <- as.vector(pGCexpdf$Included)
Totses <- as.vector(pGCexpdf$Total.Session)
pGCexpdf$freq <- as.vector(Inc/Totses)
tem26<- pGCexpdf$freq 
pGCexpdf$l10freq <- as.vector(sapply(tem26, function (x) {log(x,base = 10)}))

gPM3<-get_parameters(pGCbinomial.M2)
paramver.M2 <-as.vector(gPM3[,2]) 


df.fitted.meansM2 <- as.data.frame(matrix(data=NA, nrow=1, ncol=18))
df.fitted.meansM2[1,] <- fitted(pGCbinomial.M2)
df.fitted.meansM2
df.fitted.meansM2<- rbind(df.fitted.meansM2, rep(1,18))
df.fitted.meansM2<- rbind(df.fitted.meansM2, c(pGCexpdf$Total.Session))
df.fitted.meansM2<- rbind(df.fitted.meansM2,(fitted(pGCbinomial.M2)/pGCexpdf$Total.Session))
df.fitted.meansM2
rownames(df.fitted.meansM2)<-c("FIT","grp","Total.Session","fitted.total")
colnames(df.fitted.meansM2)<-c("N1","N2","N3","N4","C1","C2","C3","C4","N5","N6","N7","C5","C6","C7","C8","C9","C10","C11")
df.fitted.meansM2<-as.data.frame(t(df.fitted.meansM2))
df.fitted.meansM2
###

###dataframe with parameter details


df.effects.M2 <- as.data.frame(matrix(data=NA, nrow=4, ncol=13), rownames = c("A","B") )
df.effects.M2[,1]  <- c("20C","10C","20C","10C")
df.effects.M2[,2]  <- c(1.907217,4.694687,2.751394,7.180109)  #
df.effects.M2[,3]  <- c(0.9092426,3.028811,1.5237252,5.525665) #
df.effects.M2[,4]  <- c(4.000555,7.276811,4.968200,9.329912) #
df.effects.M2[,5] <-  rep(2.725086,4) #offset
df.effects.M2[,6] <- sapply(df.effects.M2[,5], function(x) {log(x,base=10)})
df.effects.M2[1,7] <- 7/56
df.effects.M2[2,7] <- 20/65
df.effects.M2[3,7] <- 11/61
df.effects.M2[4,7] <- 56/119
df.effects.M2[,8] <- sapply(df.effects.M2[,7], function(x) {log(x,base=10)})
df.effects.M2[,9] <- c(1,1,1,1)
df.effects.M2[1,10] <- df.fitted.meansM2[[1,"fitted.total"]]
df.effects.M2[2,10] <- df.fitted.meansM2[[5,"fitted.total"]]
df.effects.M2[3,10] <- df.fitted.meansM2[[9,"fitted.total"]]
df.effects.M2[4,10] <- df.fitted.meansM2[[12,"fitted.total"]]
df.effects.M2[1,11] <- exp(paramver.M2[1])
df.effects.M2[2,11] <- exp(paramver.M2[1] + paramver.M2[2])
df.effects.M2[3,11] <- exp(paramver.M2[1] + paramver.M2[3])
df.effects.M2[4,11] <- exp(paramver.M2[1] + paramver.M2[2] + paramver.M2[3])
df.effects.M2[,12] <- c("EPM","EPM","GPT","GPT")
df.effects.M2[,13] <- c("1","2","3","4")
colnames(df.effects.M2) <- c("Temperature","effects","lower","upper","offset","log(effect)","wmean", "lwmean","grp","fitted.total","fitted.from.param","Experiment","Group")

############## put both dataframes in numeric and character
pGCexpdf[,1]<-as.character(pGCexpdf[,1])
pGCexpdf[,2]<-as.character(pGCexpdf[,2])
pGCexpdf[,3]<-as.numeric(pGCexpdf[,3])
pGCexpdf[,4]<-as.numeric(pGCexpdf[,4])
pGCexpdf[,5]<-as.character(pGCexpdf[,5])
pGCexpdf[,6]<-as.character(pGCexpdf[,6])
pGCexpdf[,7]<-as.numeric(pGCexpdf[,7])
pGCexpdf[,8]<-as.numeric(pGCexpdf[,8])
##cols 9 to 11 are numeric already

####### fusion of dataframes in one to have just one for the whole plot
df.effects.M2
pGCexpdf


confintM2<- as.data.frame(confint.default(pGCbinomial.M2)) ## puts the confidence intervals of effect values to a dataframe. The effects (and therefore, the intervals) are relative to the intercept (greenhouse treatment)
confintM2adj<- data.frame(matrix(data=NA,nrow=3,ncol=2)) 
confintM2adj[1,] <- confintM2[1,]
confintM2adj[2,] <- (confintM2[1,] + confintM2[2,])
confintM2adj[3,] <- (confintM2[1,] + confintM2[3,])
## this is the transformation of the CIs of the effects (relative to greenhouse) to CIs of the frequencies (fixed), expressed still in natural log
confintM2adj
confintM2lower <-as.vector(confintM2adj[,1]) ## this puts both CIs ordered their own dataframe
confintM2upper <-as.vector(confintM2adj[,2]) ##

cisat<-add_ci(pGCexpdf,pGCbinomial.M2,alpha=0.05)
dataPGCplusCI<-as.data.frame(cisat)
dataPGCplusCI$Obsrate <- dataPGCplusCI[,"Included"]/dataPGCplusCI[,"Total.Session"]
dataPGCplusCI$Fittedrate <- fitted(pGCbinomial.M2)/dataPGCplusCI[,"Total.Session"]
dataPGCplusCI$logObsrate <- sapply(dataPGCplusCI["Obsrate"], function(x) {log(x, base=10)})
dataPGCplusCI$logTotalSeeds <- sapply(dataPGCplusCI["Total.Session"], function(x) {log(x, base=10)})
dataPGCplusCI$logFittedrate <- sapply(dataPGCplusCI["Fittedrate"], function(x) {log(x,base=10)})
dataPGCplusCI$Group <- as.character(dataPGCplusCI[,"Group"])
dataPGCplusCI  

df.for.plotting.M2 <- as.data.frame(bind_rows(dataPGCplusCI, df.effects.M2))  # Bind as rows
df.for.plotting.M2

mycolorsp7 <- c("#ebb447", "#5281e0")
names(mycolorsp7) <- c("20C","10C")

plotM2 <- ggplot(df.for.plotting.M2, aes(x=Experiment, y= Obsrate,color= Temperature)) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_jitter(aes(size = Total.Session), position = position_dodge2(width= 0.2) , alpha=0.5, show.legend= TRUE) + 
  geom_ribbon(aes( ymin=LCB0.025, ymax=UCB0.975), linetype=1, show.legend= FALSE, size=1) +  #plots the CI
  geom_errorbar(aes(x= Experiment, ymin = pred, ymax = pred), col="#000000", linetype=1, size=0.1) + #Plots the mean per group
  scale_color_manual(values= mycolorsp7) +
  scale_fill_manual(values= mycolorsp7) +
  theme(axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(0,1), expand = c(0,0))
plotM2


###################################### Model 3 #####################################################
#Necessary file: data_for_model_3.txt

#ExperimentalGroup Subject SubgColdWT SubgGHWT Experiment Temperature Genotype Count Total.seeds
#1 D z R1 A A A 16 942992 
#1 E z R2 A A A 6 238423
#1 F z R3 A A A 10 637379
#1 G z R4 A A A 1 275794
#4 H R1 z A B A 286 167169
#4 I R2 z A B A 402 147661
#4 J R3 z A B A 463 164400
#2 K z z B A A 0 14060 
#3 L z z B A B 24 18060
#5 M z z B B A 16 9560
#7 N z z B B B 198 9060 
#6 O z z C B A 31 12000
#8 P z z C B B 396 12500

### in Experiment column, A= Experiment 1 , B= Experiment 2, and C= Experiment 3
### in Genotype column, A is wt, B is dpd1 mutant
### in Temperature column, A is greenhouse condition (20-25C) , and B is chilling stress (10)
### z are empty cells

cold.dpd1.exp1 <- read.table("data_for_model_3.txt", header = TRUE, na.strings="z")
TGRexp <- as.data.frame(cold.dpd1.exp1)


cold.dpd1.nb.M3 <- glm.nb(Count ~ Genotype + Temperature + Experiment + Genotype:Temperature + offset(log(Total.seeds)), contrasts=list(Genotype=contr.treatment,
                                                                                                                                        Temperature=contr.treatment, 
                                                                                                                                        Experiment=contr.treatment(3,base=1)), data=TGRexp)

cold.dpd1.poisson.m3.2 <- glm(Count ~ Genotype + Temperature + Experiment + Genotype:Temperature + offset(log(Total.seeds)), family= poisson(link="log"), 
                              contrasts=list(Genotype=contr.treatment,
                                             Temperature=contr.treatment,                                                                                                            
                                             Experiment=contr.treatment(3,base=1)), data=TGRexp) 

cold.dpd1.nb.m3.3 <- glm.nb(Count ~ Genotype + Temperature + Experiment + Genotype:Temperature + Genotype:Experiment + offset(log(Total.seeds)), contrasts=list(Genotype=contr.treatment,
                                                                                                                                                                Temperature=contr.treatment, 
                                                                                                                                                                Experiment=contr.treatment(3,base=1)), data=TGRexp)

cold.dpd1.nb.m3.4 <- glm.nb(Count ~ Genotype + Temperature + Experiment + Genotype:Temperature + Genotype:Experiment + Temperature:Experiment + offset(log(Total.seeds)), contrasts=list(Genotype=contr.treatment,
                                                                                                                                                                Temperature=contr.treatment, 
                                                                                                                                                                Experiment=contr.treatment(3,base=1)), data=TGRexp)

df.mastertable.M3<-as.data.frame(add_ci(TGRexp,cold.dpd1.nb.M3,alpha=0.05))
df.mastertable.M3[14,] <- c(9,NA,NA,NA,"B","B","B",1,100000,NA,NA,NA) #14th value is just for point scale purposes
df.mastertable.M3$Count <-as.numeric(df.mastertable.M3[,"Count"])
df.mastertable.M3$Total.seeds <-as.numeric(df.mastertable.M3[,"Total.seeds"])
df.mastertable.M3$pred <-as.numeric(df.mastertable.M3[,"pred"])
df.mastertable.M3$LCB0.025 <-as.numeric(df.mastertable.M3[,"LCB0.025"])
df.mastertable.M3$UCB0.975 <-as.numeric(df.mastertable.M3[,"UCB0.975"])
df.mastertable.M3$Obsrate <- df.mastertable.M3[,"Count"]/df.mastertable.M3[,"Total.seeds"]
df.mastertable.M3$Predrate <- df.mastertable.M3[,"pred"]/df.mastertable.M3[,"Total.seeds"]
df.mastertable.M3$Fittedrate <- fitted(cold.dpd1.nb.M3)/df.mastertable.M3[,"Total.seeds"]
df.mastertable.M3[14,"Fittedrate"] <- NA
df.mastertable.M3$LCBrate <- df.mastertable.M3[,"LCB0.025"]/df.mastertable.M3[,"Total.seeds"]
df.mastertable.M3$UCBrate <- df.mastertable.M3[,"UCB0.975"]/df.mastertable.M3[,"Total.seeds"]
df.mastertable.M3$logObsrate <- sapply(df.mastertable.M3["Obsrate"], function(x) {log(x, base=10)})
df.mastertable.M3$logTotalSeeds <- sapply(df.mastertable.M3["Total.seeds"], function(x) {log(x, base=10)})
df.mastertable.M3$logPredrate <- sapply(df.mastertable.M3["Predrate"], function(x) {log(x,base=10)})
df.mastertable.M3$logFittedrate <- sapply(df.mastertable.M3["Fittedrate"], function(x) {log(x,base=10)})
df.mastertable.M3$logLCBrate <- sapply(df.mastertable.M3["LCBrate"], function(x) {log(x,base=10)})
df.mastertable.M3$logUCBrate <- sapply(df.mastertable.M3["UCBrate"], function(x) {log(x,base=10)})
df.mastertable.M3$Condition <- c(1,1,1,1,2,2,2,1,3,2,4,2,4,NA)#14th value is just for point scale purposes
df.mastertable.M3$Condition2<- as.character(df.mastertable.M3[,"Condition"])
df.mastertable.M3$ExperimentalGroup <- as.character(df.mastertable.M3$ExperimentalGroup)
df.mastertable.M3[8,"logObsrate"] <- -7 ### this is only so the point with 0 obs appears in the plot, is going to be put in its own zone outside of R
df.mastertable.M3[14,"logObsrate"] <- -7 
df.mastertable.M3$WeiMeanRate <- c(-4.80549,-4.80549,-4.80549,-4.80549,-2.6212,-2.6212,-2.6212,-4.80549,-2.87651,-2.6212,-1.55986,-2.6212,-1.55986,NA)

plotM3 <- ggplot(data=df.mastertable.M3, aes(x = ExperimentalGroup, y = logObsrate, color= Condition2)) +
  geom_point(aes(size = logTotalSeeds, shape=Experiment), position = position_dodge2(width=1), alpha=1, show.legend= TRUE) +
  scale_size_continuous(range = c(4,8)) +
  scale_shape_manual(values=c(1, 10, 13)) +
  geom_ribbon(aes(x=ExperimentalGroup, ymin=logLCBrate, ymax=logUCBrate), linetype=1, show.legend= FALSE, size=1) + 
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black"))  +
  geom_errorbar(aes(x= ExperimentalGroup, ymin = logPredrate, ymax = logPredrate, color= Condition2), linetype=1, size=1.3, alpha=0.5)+ #Plots the mean per group
  geom_errorbar(aes(x= ExperimentalGroup, ymin = WeiMeanRate, ymax = WeiMeanRate), col="#000000", linetype=1, size=1) + #Plots the mean per group
  scale_x_discrete(position = "top") +
  scale_y_continuous(limits=c(-7,0), expand = c(0,0)) +
  xlab("Group") +
  ylab("Rate of PPI")
plotM3


################################################### Model 4 #####################################################################################

#Temperature Total.pollen Alive
#20C 747 628 
#20C 367 310
#20C 376 314
#10C 392 87
#10C 622 134
#10C 643 126

pollen.survival.1 <- read.table("data_for_model_4.txt", header = TRUE)
PSexp1 <- as.data.frame(pollen.survival.1)
PSexp1$Freq <- PSexp1$Alive/PSexp1$Total.pollen
PSexp1$Temperature = factor(PSexp1$Temperature, levels=unique(PSexp1$Temperature))




PSexp1bin.M4 <- glm(Freq ~ Temperature, weights = Total.pollen, family=binomial(link="log"),
                    data=PSexp1, contrasts = contr.treatment)





PSexp1$l10freq <- sapply(PSexp1$Freq, function (x) {log(x,base = 10)})
gpmM4<-get_parameters(PSexp1bin.M4)
paramver.M4 <-as.vector(gpmM4[,2]) 



df.fitted.M4<-as.data.frame(add_ci(PSexp1,PSexp1bin.M4,alpha=0.05))

df.effects.M4 <- as.data.frame(matrix(data=NA, nrow=2, ncol=12), rownames = c("A","B") )
df.effects.M4[1,1]  <- c("20C")
df.effects.M4[2,1]  <- c("10C")
df.effects.M4[,1]  <- as.factor(df.effects.M4[,1])
df.effects.M4[1,2] <-  1252/1490
df.effects.M4[2,2] <-  347/1657
df.effects.M4[1,3] <- log(1252/1490, base = 10)
df.effects.M4[2,3] <- log(347/1657, base = 10)
df.effects.M4[,4] <- c(1,2)
df.effects.M4[1,5] <- df.fitted.M4[[1,"pred"]]
df.effects.M4[2,5] <- df.fitted.M4[[4,"pred"]]
df.effects.M4[1,6] <- exp(paramver.M4[1])
df.effects.M4[2,6] <- exp(paramver.M4[1] + paramver.M4[2])
colnames(df.effects.M4) <- c("Temperature","wmean", "l10wmean","grp","fitted.total","fitted.from.param")

df.for.plotting.M4 <- as.data.frame(dplyr::bind_rows(df.fitted.M4, df.effects.M4))  # Bind as rows

mycolorsp8 <- c("#ebb447", "#5281e0")
names(mycolorsp8) <- c("20C","10C")
plotM4 <- ggplot(df.for.plotting.M4 , aes(x=factor(Temperature, level=c("20C","10C")), y=Freq, color= Temperature)) +
  scale_color_manual(values= mycolorsp8) +
  scale_fill_manual(values= mycolorsp8) +
  geom_errorbar(aes(x=factor(Temperature, level=c("20C","10C")), ymin = fitted.total, ymax = fitted.total), col="#000000", linetype=1, size=0.1) + #Plots the mean per group
  geom_ribbon(aes(x=factor(Temperature, level=c("20C","10C")), ymin=LCB0.025, ymax=UCB0.975), linetype=1, show.legend= FALSE, size=1) +
  geom_point(aes(size= Total.pollen) , position = position_dodge2(width=0.7), alpha=0.5, show.legend= TRUE) +
  scale_size_continuous(range = c(5,10))  +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
  xlab("Temperature") +
  ylab("Proportion of viable pollen")
plotM4

################################################ Model 5 ##########################################################

#Point Replicate Genotype Total.pollen Alive 
#AA C A 58 51
#AB C A 81 69
#AC C A 92 79
#AD C A 47 42
#AE C A 45 41
#AF D A 86 74
#AG D A 76 69
#AH D A 90 78
#AI D A 87 73
#AJ D A 35 26
#AK E A 42 36
#AL E A 35 26
#AM E A 36 28
#AN E A 24 19
#AO E A 24 17
#AP E A 20 17
#AQ E A 29 21
#AR E A 33 27
#AS E A 37 29
#T F B 62 45
#U F B 76 46
#V F B 81 55
#W F B 60 40
#X F B 56 40
#Y G B 104 71
#Z G B 96 56
#ZA G B 128 90
#ZB G B 183 132
#ZC G B 105 65
#ZD H B 88 43
#ZE H B 63 42
#ZF H B 48 24
#ZG H B 47 32
#ZH H B 54 35

### A and B are wt and dpd1 genotype
### Replicates C, D and E are from wt. F G and H from dpd1 

PSexpdpd1 <- read.table("data_for_model_5.txt", header = TRUE)
Polsurdpd1 <- as.data.frame(PSexpdpd1)
Polsurdpd1$Freq <- Polsurdpd1$Alive / Polsurdpd1$Total.pollen 
Polsurdpd1$Genotype = factor(Polsurdpd1$Genotype, levels=unique(Polsurdpd1$Genotype))
Polsurdpd1$Replicate = factor(Polsurdpd1$Replicate, levels=unique(Polsurdpd1$Replicate))


PSexp2.bin.M5 <- glm(Freq ~ Genotype + Replicate, weights= Total.pollen, family=binomial(link="log"),
                     contrasts=list(Genotype=contr.treatment(2,base=1),Replicate=contr.treatment(6,base=1)),data=Polsurdpd1)

PSexp2.bin.m5.2 <- glm(Freq ~ Genotype, weights= Total.pollen, family=binomial(link="log"),
                       contrasts=list(Genotype=contr.treatment(2,base=2)), data=Polsurdpd1)

Polsurdpd1$freq <- Polsurdpd1$Alive / Polsurdpd1$Total.pollen
Polsurdpd1$l10freq <- sapply(Polsurdpd1$freq, function (x) {log(x,base = 10)})

gpM5<-get_parameters(PSexp2.bin.M5) ## this function extracts the parameter values from the model
paramver.M5 <-as.vector(gpM5[,2])

df.fitted.means.M5 <- as.data.frame(matrix(data=NA, nrow=1, ncol=34))
df.fitted.means.M5[1,] <- fitted(PSexp2.bin.M5)
df.fitted.means.M5<- rbind(df.fitted.means.M5, rep(1,6))
df.fitted.means.M5<- rbind(df.fitted.means.M5, c(Polsurdpd1$Total.pollen))
rownames(df.fitted.means.M5)<-c("fitted.total","grp","Total.pollen")
colnames(df.fitted.means.M5)<-c("C","D","E","F","H","G")
df.fitted.means.M5<-as.data.frame(t(df.fitted.means.M5))


df.effects.M5 <- as.data.frame(matrix(data=NA, nrow=6, ncol=6), rownames = c("A","B") )
df.effects.M5[,1]  <- c("A","A","A","B","B","B")
df.effects.M5[,2]  <- c("C","D","E","F","G","H")
df.effects.M5[,3] <-  c(822/977,822/977,822/977,816/1251,816/1251,816/1251)
df.effects.M5[,4] <- c(rep(log(822/977, base = 10),3), rep(log(816/1251, base = 10),3))
df.effects.M5[,5] <- c(1,2,3,4,5,6)
df.effects.M5[1,6] <- df.fitted.means.M5[[1,"fitted.total"]]
df.effects.M5[2,6] <- df.fitted.means.M5[[6,"fitted.total"]]
df.effects.M5[3,6] <- df.fitted.means.M5[[11,"fitted.total"]]
df.effects.M5[4,6] <- df.fitted.means.M5[[20,"fitted.total"]]
df.effects.M5[5,6] <- df.fitted.means.M5[[25,"fitted.total"]]
df.effects.M5[6,6] <- df.fitted.means.M5[[30,"fitted.total"]]
colnames(df.effects.M5) <- c("Genotype","Replicate","wmean", "l10wmean","grp","fitted.total")


df.for.plotting.M5 <- as.data.frame(dplyr::bind_rows(Polsurdpd1, df.effects.M5))  # Bind as rows

mycolorsp8 <- c("#ebb447", "#c77cff")
names(mycolorsp8) <- c("A","B")
plotM5 <- ggplot(df.for.plotting.M5, aes(x=Replicate, y=freq, color= Genotype)) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(aes(size= Total.pollen), alpha=0.5,position = position_dodge2(width=0.55), show.legend= TRUE) +
  scale_size_continuous(range = c(4,7))  +
  geom_errorbar(aes(x = Replicate, ymin = wmean, ymax = wmean), col="#000000", linetype=1, size=1) +
  geom_errorbar(aes(x = Replicate, ymin = fitted.total, ymax = fitted.total), col="#000000", linetype=1, size=0.2) +
  scale_color_manual(values= mycolorsp8) +
  scale_fill_manual(values= mycolorsp8) +
  theme(axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(0,1), expand = c(0,0))+
  xlab("Genotype and Replicate") +
  ylab("Proportion of viable pollen")
plotM5

########################################################Models, Parameters (untransformed), AICc, log Likelihood#########################

sink('Models_untransformed_parameters_AICc_logLik.txt')

########## Parameters and SEs in these reports are in log e, they have to be transformed to match those in Extended Data Table 2

print("M1 negative binomial: LS.exp1.nb.M1") 
summary(LSexp1.nb.M1)
print("AICc:")
AICc(LSexp1.nb.M1)
logLik(LSexp1.nb.M1)

print("m1.2 poisson: LS.exp1.poisson.m1.2")

summary(LSexp1.poisson.m1.2)
logLik(LSexp1.poisson.m1.2)
print("AICc:")
AICc(LSexp1.poisson.m1.2)

### 

print("M2 binomial (without interaction term): pGCbinomial.M2")

summary(pGCbinomial.M2)
print("AICc:")
AICc(pGCbinomial.M2)
logLik(pGCbinomial.M2)

print("m2.2 binomial (with interaction term): pGCbinomial.m2.2") 

summary(pGCbinomial.m2.2)
print("AICc:")
AICc(pGCbinomial.m2.2)
logLik(pGCbinomial.m2.2)

###

print("M3 Negative binomial: cold.dpd1.nb.M3")

summary(cold.dpd1.nb.M3)
print("AICc:")
AICc(cold.dpd1.nb.M3)
logLik(cold.dpd1.nb.M3)

print("m3.2 Poisson: cold.dpd1.nb.m3.2")

summary(cold.dpd1.poisson.m3.2)
print("AICc:")
AICc(cold.dpd1.poisson.m3.2)
logLik(cold.dpd1.poisson.m3.2)

print("m3.3 Negative binomial (with GxE2 interaction term):  cold.dpd1.nb.m3.3")

summary(cold.dpd1.nb.m3.3)
print("AICc:")
AICc(cold.dpd1.nb.m3.3)
logLik(cold.dpd1.nb.m3.3)

print("m3.4 Negative binomial (with GXE2 and TxE interaction terms):  cold.dpd1.nb.m3.4")

summary(cold.dpd1.nb.m3.4)
print("AICc:")
AICc(cold.dpd1.nb.m3.4)
logLik(cold.dpd1.nb.m3.4)

print("M4 binomial: PSexp1bin.M4")

summary(PSexp1bin.M4)
print("AICc:")
AICc(PSexp1bin.M4)
logLik(PSexp1bin.M4)

print("M5 binomial: PSexp2.bin.M5")

summary(PSexp2.bin.M5)
print("AICc:")
AICc(PSexp2.bin.M5)
logLik(PSexp2.bin.M5)

print("m5.2 binomial (without replicate parameters): PSexp2.bin.m5.2")

summary(PSexp2.bin.m5.2)
print("AICc:")
AICc(PSexp2.bin.m5.2)
logLik(PSexp2.bin.m5.2)

sink()

###################### Likelihood Ratio Tests ########################

sink('Likelihood_ratio_test.txt')

### Residual deviance values can be found in the summary of the models, directly above,
### but the degrees of freedom shown of said deviance on are incorrect in the case of negative binomial models. 
### -1df should be substracted on those. The error arises because the output does not count the overdispersion parameter in the calculation 
### Correct degrees of freedom are shown in the LRT part of Extended Data Table 1 

####### LRT overdispersion tests 

print("LRT overdispersion test")
print("M1 vs m1.2")
odTest(LSexp1.nb.M1) 


print("M3 vs m3.2")
odTest(cold.dpd1.nb.M3)

######## Likelihood ratio test of saturated models

print("M1 vs saturated M1")
print("residual deviance: 16.65, df: 10")
print( "1 - pchisq(16.65,10)")
1 - pchisq(16.65,10)

print("M2 vs saturated M2")
print("residual deviance: 18.63, df: 15")
print( "1 - pchisq(18.63,15)")

1 - pchisq(18.63,15)

print("M3 vs saturated M3")
print("residual deviance: 11.27, df: 6")
print( "1 - pchisq(11.27,6)")

1 - pchisq(11.27,6)

print("M4 vs saturated M4")
print("residual deviance: 1.345, df: 4")
print( "1 - pchisq(1.345,4)")

1 - pchisq(1.345,4)


print("M5 vs saturated M5")
print("residual deviance: 29.83, df: 4")
print( "1 - pchisq(29.83,28)")

1 - pchisq(29.83,28)

######## Comparison between M5 and M5.2

print("M5 residual deviance: 29.83, df= 28")
print("m5.2 residual deviance: 46.34, df= 32")
print("Likelihood ratio: 46.34 - 29.83, deltadf=4")
print("Likelihood ratio 16.51, deltadf=4")
print( "1 - pchisq(16.51,4)")

1 - pchisq(16.51,4)

sink()

#################### Reports of Parameter values, Standard Error, Wald CI95 (transformed and untransformed) z-value, p-values #########################

main.model.list <- as.vector(c(LSexp1.nb.M1,pGCbinomial.M2,cold.dpd1.nb.M3,PSexp1bin.M4,PSexp2.bin.M5))
SM1<-summary(LSexp1.nb.M1)
SM2<-summary(pGCbinomial.M2)
SM3<-summary(cold.dpd1.nb.M3)
SM4<-summary(PSexp1bin.M4)
SM5<-summary(PSexp2.bin.M5)

PM1 <- get_parameters(LSexp1.nb.M1)
PM2 <- get_parameters(pGCbinomial.M2)
PM3 <- get_parameters(cold.dpd1.nb.M3)
PM4 <- get_parameters(PSexp1bin.M4)
PM5 <- get_parameters(PSexp2.bin.M5)
PM5 <- PM5 %>%  na.omit()

CIM1 <- confint.default(LSexp1.nb.M1)
CIM2 <- confint.default(pGCbinomial.M2)
CIM3 <- confint.default(cold.dpd1.nb.M3)
CIM4 <- confint.default(PSexp1bin.M4)
CIM5 <- confint.default(PSexp2.bin.M5)
CIM5 <- CIM5 %>%  na.omit()

PM1[,3]<-SM1$coefficients[,2]
PM2[,3]<-SM2$coefficients[,2]
PM3[,3]<-SM3$coefficients[,2]
PM4[,3]<-SM4$coefficients[,2]
PM5[,3]<-SM5$coefficients[,2]

PM1 <- cbind(PM1,CIM1)
PM2 <- cbind(PM2,CIM2)
PM3 <- cbind(PM3,CIM3)
PM4 <- cbind(PM4,CIM4)
PM5 <- cbind(PM5,CIM5)

PM1[,6]<-SM1$coefficients[,3]
PM2[,6]<-SM2$coefficients[,3]
PM3[,6]<-SM3$coefficients[,3]
PM4[,6]<-SM4$coefficients[,3]
PM5[,6]<-SM5$coefficients[,3]
PM1[,7]<-SM1$coefficients[,4]
PM2[,7]<-SM2$coefficients[,4]
PM3[,7]<-SM3$coefficients[,4]
PM4[,7]<-SM4$coefficients[,4]
PM5[,7]<-SM5$coefficients[,4]
Master.table.PM.lognatural<-as.data.frame(bind_rows(PM1,PM2,PM3,PM4,PM5))
colnames(Master.table.PM.lognatural) <- c("Parameter","Estimate","Standard.Error","CI95.lower.bound","CI95.upper.bound", "Wald.z-statistic", "P-value")
rownames(Master.table.PM.lognatural) <- c("M1 Intercept","M1 High Light", "M1 Heat", "M1 Drought", "M1 Cold", 
                                          "M2 Intercept", "M2 Cold", "M2 Experiment GPT", 
                                          "M3 Intercept", "M3 Genotype dpd1", "M3 Cold", "M3 Experiment 2", "M3 Experiment 3", "M3 Interaction dpd1/cold",
                                          "M4 Intercept", "M4 cold",
                                          "M5 Intercept", "M5 Genotype", "M5 Replicate 2","M5 Replicate 3", "M5 Replicate 5", "M5 Replicate 6")
Master.table.PM.lognatural

Master.table.PM.log10 <-Master.table.PM.lognatural
Master.table.PM.log10[,2] <-0.4342944*Master.table.PM.lognatural[,2]
Master.table.PM.log10[,3] <-0.4342944*Master.table.PM.lognatural[,3]
Master.table.PM.log10[,4] <-0.4342944*Master.table.PM.lognatural[,4]
Master.table.PM.log10[,5] <-0.4342944*Master.table.PM.lognatural[,5]
Master.table.PM.log10

sink(file= "Parameter tables log natural.txt") ### All values in base e (without transformations)
Master.table.PM.lognatural
sink()

sink(file= "Parameter tables log 10.txt") ##### Prints all the values in log 10, as listed in Extended Data Table 2
Master.table.PM.log10
sink()

######### Plots to PDF ###############

pdf(file="Model 1.pdf", height = 5, width =  8)
plotM1 #here goes the plot code###
dev.off()

pdf(file="Model 2.pdf", height = 5, width =  6)
plotM2 #here goes the plot code###
dev.off()

pdf(file="Model 3.pdf", height = 5, width =  8)
plotM3 #here goes the plot code###
dev.off()

pdf(file="Model 4.pdf", height = 5, width =  5)
plotM4 #here goes the plot code###
dev.off()

pdf(file="Model 5.pdf", height = 5, width =  7)
plotM5 #here goes the plot code###
dev.off()

############################################### END #########################################################################################################
