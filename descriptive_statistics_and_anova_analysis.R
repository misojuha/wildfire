library(devtools)
library(ggbiplot)
library(tidyverse)
library(RNHANES)
library(dplyr)
library(MASS)

data <- read.table("data.csv", sep =",", header=T)

#rounding to 4 decimal places
data<-data %>% 
  mutate_if(is.numeric, round, digits=4)


#BOXPLOTS
boxplot(data$OM)


#QQ PLOTS
qqnorm(data$OM, pch = 1, frame = FALSE, main="")
qqline(data$OM, col = "red", lwd = 2)
qqnorm(data$pH, pch = 1, frame = FALSE, main="")
qqline(data$pH, col = "red", lwd = 2)
qqnorm(data$EC, pch = 1, frame = FALSE, main="")
qqline(data$EC, col = "red", lwd = 2)
qqnorm(data$CaCO3, pch = 1, frame = FALSE, main="")
qqline(data$CaCO3, col = "red", lwd = 2)
#Histograms
hist(data$OM, xlab="soil OM (%)", main="")
hist(data$EC,xlab="EC (uS/cm)", main="")
hist(data$CaCO3,xlab="CaCO3 (%)", main="")
hist(data$pH,xlab="pH (-log[H+])", main="")

#Shapiro wilk test for normality if p<0.05 data is not normal
shapiro.test(data$OM)
shapiro.test(data$EC)
shapiro.test(data$pH)
shapiro.test(data$CaCO3)

#BOXCOX TRANSFORMATION and q-q plot,histogram for normalized data
#OM
boxcox(lm(data$OM ~ 1))
BC_OM <- (data$OM ^ -0.948031 - 1) / -0.948031
qqnorm(BC_OM, pch = 1, frame = FALSE, main="")
qqline(BC_OM, col = "red", lwd = 2)
hist(BC_OM,xlab="soil OM(%)_transformed data", main="")
data<-bind_cols(data,BC_OM)
#renaming the new column to BC_OM ...12 was old column name
data<-rename(data, BC_OM = ...12)
#remove 1 column
data <- subset (data, select = -...13)
#EC
BC_EC <- (data$EC ^ -0.673821 - 1) / -0.673821
qqnorm(BC_EC, pch = 1, frame = FALSE, main="")
qqline(BC_EC, col = "red", lwd = 2)
hist(BC_EC,xlab="EC (??S/cm)_transformed data", main="")

#Yeo-Johnson transformation for CaCO3 (which contains zero values)
qqnorm(data$CaCO3, pch = 1, frame = FALSE, main="")
qqline(data$CaCO3, col = "red", lwd = 2)
install.packages("VGAM")
library(VGAM)
lambda.max <- -10.9598 #the optimal lambda value
yj<-yeo.johnson(data$CaCO3, lambda.max, derivative = 0,
            epsilon = sqrt(.Machine$double.eps), inverse = FALSE)
qqnorm(yj, pch = 1, frame = FALSE, main="")
qqline(yj, col = "red", lwd = 2)
yj
#alternative Y-J transformation of CaCO3 (result is the same) but this uses car package
CaCO3_t<-yjPower(data$CaCO3, lambda=lambda.max, jacobian.adjusted=FALSE)
qqnorm(CaCO3_t, pch = 1, frame = FALSE, main="")
qqline(CaCO3_t, col = "red", lwd = 2)


#CALCULATE Z-SCORES FOR OUTLIER DETECTION: cutoff at +-3 std
z_scores_OM <- (data$OM-mean(data$OM))/sd(data$OM)
x_zscore <- (data$OM - mean(data$OM, na.rm = TRUE)) / sd(data$OM, na.rm = TRUE)
x_zscore
z_scores_OM


#outlier removal
data2 <- data[-c(8), ]


#2 way ANOVA
anova.OM2way<-lm(OM~Time+Severity+Time:Severity,data=data)
summary(anova.OM2way)
data2 <- data[-c(485,495,521,525,538), ]
shapiro.test(data2$OM)
qqnorm(data2$OM, pch = 1, frame = FALSE)
qqline(data2$OM, col = "steelblue", lwd = 2)
anova.OM2way2<-lm(OM~Time+Severity+Time:Severity,data=data)
summary(anova.OM2way2)

anova.pH2way<-lm(pH~Time+Severity+Time:Severity,data=data)

install.packages("performance")
install.packages("see")
install.packages("patchwork")
library("performance")
library("see")
library("patchwork")
install.packages("effects")
install.packages("sjPlot")
library("sjPlot")


#OM temporal distribution

plot_model(anova.OM2way, type = "pred", terms = c("Severity", "Time"))
plot_model(anova.OM2way, type = "pred", terms = c("Time", "Severity"),colors = c("dark green", "orange", "red"))
plot_model(anova.OM2way2, type = "pred", terms = c("Time", "Severity"),colors = c("Dark2"))
#pH temporal distribution
plot_model(anova.pH2way, type = "pred", terms = c("Time", "Severity"),colors = c("dark green", "orange", "red"))
#LETTERS for OM~time
ANOVAOM=aov(anova.OM2way2)

# Tukey test to study each pair of treatment 
TUKEY <- TukeyHSD(aov(OM ~ Severity+Time+Severity:Time, data = data2))

# Tuckey test representation :
plot(TUKEY , las=1 , col="brown")

library(lsmeans)

marginal = lsmeans(ANOVAOM,
                   ~ Time)

pairs(marginal,
      adjust="tukey")
CLD = multcomp::cld(marginal,
                    alpha   = 0.05,
                    Letters = letters,    ### Use lower-case letters for .group
                    adjust  = "tukey")
CLD

#LETTERS for EC~Time
anova.EC2way<-lm(EC~Time+Severity+Time:Severity,data=data2)
ANOVAEC=aov(anova.EC2way)

# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(aov((EC^-0.675262-1/-0.675262) ~ Severity+Vegetation+Severity:Vegetation, data = data2))


marginal = lsmeans(ANOVAEC,
                   ~ Time)

pairs(marginal,
      adjust="tukey")
CLD = multcomp::cld(marginal,
                    alpha   = 0.05,
                    Letters = letters,    ### Use lower-case letters for .group
                    adjust  = "tukey")
CLD



# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'Severity', conf.level=0.95)
TUKEY

# Tukey test representation :
plot(TUKEY , las=1 , col="brown")

#___________________
#___________________
#___________________

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)


#________________________________
#________________________________
#________________________________

#3 way anova
#OM
three.way <- aov(OM ~ Time + Severity+Vegetation+Time:Severity+Time:Vegetation+Severity:Vegetation, data = data2)
summary(three.way)
TukeyHSD(aov(OM ~ Severity+Severity:Vegetation, data = data))

#nice boxplot for use on congress presentation
anova.OM<-lm(OM ~ Time + Severity+Vegetation+Time:Severity+Time:Vegetation+Severity:Vegetation,data=data)
plot_model(anova.OM)
plot_model(anova.OM, type = "pred", terms = c("Severity", "Time", "Vegetation"))
plot_model(anova.OM, type = "pred", terms = c("Time", "Severity", "Vegetation"), colors = "Dark2")

plot_model(anova.OM, type = "pred", terms = c("Time", "Severity", "Vegetation"), colors =c("dark green", "orange", "red"))
#pH, EC and CaCO3
anova.pH<-lm(pH ~ Time + Severity+Vegetation+Time:Severity+Time:Vegetation+Severity:Vegetation,data=data)
plot_model(anova.pH, type = "pred", terms = c("Time", "Severity", "Vegetation"), colors =c("dark green", "orange", "red"))

anova.EC<-lm(EC ~ Time + Severity+Vegetation+Time:Severity+Time:Vegetation+Severity:Vegetation,data=data)
plot_model(anova.EC, type = "pred", terms = c("Time", "Severity", "Vegetation"), colors =c("dark green", "orange", "red"))

anova.CaCO3<-lm(CaCO3 ~ Time + Severity+Vegetation+Time:Severity+Time:Vegetation+Severity:Vegetation,data=data)
plot_model(anova.CaCO3, type = "pred", terms = c("Time", "Severity", "Vegetation"), colors =c("dark green", "orange", "red"))
###CORRELATION plot#####
#___________________________#
install.packages('corrplot')
library(corrplot)
data3 <- data2[,c(5,6,7,9)]

corrplot<-cor(data3,use="pairwise.complete.obs")

corrplot(corrplot, order='AOE')
#calculation of p level
testRes = cor.mtest(data2, conf.level = 0.95)
#add p value to corrplot
corrplot(corrplot, p.mat = testRes$p, insig = 'p-value', sig.level = -1)
corrplot(corrplot, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')


