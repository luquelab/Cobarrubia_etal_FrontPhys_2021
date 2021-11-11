sink("ResearchCode.txt")

# Preamble ----------------------------------------------------------------
#
#PURPOSE: The investigation of nanoparticles diffusivity and transport through mucus 
#or bio-like viscoelastic fluids. The anallysis of nanoparticle diffusivity is to 
#any nanoparticle characteristics that viruses, mainly bacteriophages, need to 
#possess in order to explain its mobile behavior.
#
#OBJECTIVE: This code is to  statistically analyize collected data from scientific 
#literature. The code also plots statisical correlation with diameter, charge, effective
#diffusion and subdiffusional analysis. 
#
#Version: 
#DATE       Name                        Version             Details
#11/27/18    A. Cobarrubia                  1                Creation of Code
#12/3/18    A. Cobarrubia                  2                Fixing graphical parameters
#12/8/18    A. Cobarrubia                  3                Statisical Analysis on Effective Dif as a function of alpha
#12/10/18   A. Cobarrubia                  4                Categorizing Based on Adhesive and Inert models
#5/24/2019  A. Cobarrubia                  5                Added Prediction alpha values with box plots about size,zeta,Deff,and alpha;
#                                                           Added Temp. and pH graphs. 
#6/18/19    A. Cobarrubia                  7                Added 3d plots, dominate mucin analysis and random forest analysis
#
#9/30/19  Antonio Cobarrubia               8                Rename some of the outline and variables


# Packages ----------------------------------------------------------------

library(plyr)
library(caTools)
library(combinat)
library(cowplot)
library(dplyr)
library(egg)
library(ggplotify)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(ggsignif)
library(gridExtra)
library(gtable)
library(knitr)
library(magrittr)
library(missForest)
library(oceanmap)
library(plotly)
library(plot3D)
library(pillar)

library(randomForest)
library(rattle)
library(RColorBrewer)
library(rfPermute)
library(reshape2)
library(rpart)
library(rpart.plot)
library(scales)
library(stringr)
library(tibble)
library(tidyr)

# Read_file ---------------------------------------------------------------


# Read Data File
Core_Mucus <- read.csv(file= "common_core8.csv", header = TRUE)


# Data_frames_based_on_variables ------------------------------------------

# Fitler variables based on conditions
Size1_ind = which(Core_Mucus$Diameter<100)
Size2_ind = which(c(Core_Mucus$Diameter>100, Core_Mucus$Diameter == 100))
neg_ind = which(Core_Mucus$Zeta<0)
pos_ind = which(Core_Mucus$Zeta>0)
alpha_ind = which(Core_Mucus$alpha > 0)
cooh_ind = which(Core_Mucus$Surface_Chemistry =="COOH")
amine_ind = which(Core_Mucus$Surface_Chemistry =="AMINE")
peg_ind = which(Core_Mucus$Surface_Chemistry =="PEG")
chit_ind = which(Core_Mucus$Surface_Chemistry =="CHITOSAN")
antibody_ind = which(Core_Mucus$Surface_Chemistry =="Antibody")
virus_ind = which(Core_Mucus$Surface_Chemistry == "Virus")
mucoad_ind = which(Core_Mucus$Surface_Chemistry == "COOH" | Core_Mucus$Surface_Chemistry == 'CHITOSAN')
mucoin_ind = which(Core_Mucus$Surface_Chemistry == "PEG" | Core_Mucus$Surface_Chemistry == "AMINE" )
VA_ind = which(Core_Mucus$Surface_Chemistry =="Antibody" | Core_Mucus$Surface_Chemistry == "Virus")
nonalpha_ind = which(is.na(Core_Mucus$alpha))
###
COOH <- Core_Mucus[cooh_ind, ]
AMINE <- Core_Mucus[amine_ind, ]
PEG <-  Core_Mucus[peg_ind, ]
CHITOSAN <- Core_Mucus[chit_ind, ]
ANTIBODIES <- Core_Mucus[antibody_ind, ]
Virus <- Core_Mucus[virus_ind,]
NONALPHA <- Core_Mucus[nonalpha_ind,]
Negative <- Core_Mucus[neg_ind, ]
Positive <- Core_Mucus[pos_ind, ]
Mucoadhesive <- Core_Mucus[mucoad_ind, ]
Mucoinert <- Core_Mucus[mucoin_ind, ]
Bioparticles <- Core_Mucus[VA_ind,]
alphapoints <- Core_Mucus[alpha_ind,]
SmallMucus <- Core_Mucus[Size1_ind,]
LargeMucus <- Core_Mucus[Size2_ind,]



Tissue1 <- data.frame(Mucus_Type = Core_Mucus$Mucus_Type)

TTISCM_ind <- which(Tissue1$Mucus_Type == "Human cervical mucus"                | Tissue1$Mucus_Type == "Cervicovaginal mucus" | Tissue1$Mucus_Type == "Cervicovaginal mucus " )
TTISIN_ind <- which(Tissue1$Mucus_Type == "Pig porcine ileum intestinal mucus"  | Tissue1$Mucus_Type == "Porcine intestinal mucus")                                
TTISST_ind <- which(Tissue1$Mucus_Type == "Porcine gastric mucins"              | Tissue1$Mucus_Type == "Porcine gastric mucus type II")                         
TTISCF_ind <- which(Tissue1$Mucus_Type == "Cystic fibrosis")                                                                                                         
TTISMA_ind <- which(Tissue1$Mucus_Type == "Matrigel" | Tissue1$Mucus_Type == "Matrigel-biotinylate")                                                              
TTISRS_ind <- which(Tissue1$Mucus_Type == "Respiratory mucus")                                                                                                       

TTISHL_ind <- which(Tissue1$Mucus_Type == "Respiratory mucus" | Tissue1$Mucus_Type == "Cystic fibrosis")

Tissue1[TTISCM_ind,"Mucus_Source"] <- "Human cervix"
Tissue1[TTISIN_ind,"Mucus_Source"] <- "Pig intestine"
Tissue1[TTISST_ind,"Mucus_Source"] <- "Pig stomach"
Tissue1[TTISRS_ind,"Mucus_Source"] <- "Human lung"
Tissue1[TTISCF_ind,"Mucus_Source"] <- "Human lung"
Tissue1[TTISMA_ind, "Mucus_Source"] <- "Hydrogel"
Tissue1$Mucus_Type <- NULL
Core_Mucus <- cbind(Core_Mucus,Tissue1)

Core_Mucus$Mucus_Source <- as.factor(Core_Mucus$Mucus_Source)

# Old Analysis ------------------------------------------------------------


### Anomalous Diffusion

Adiameter = lm(alpha ~ log10(Diameter), data = Core_Mucus)
summary(Adiameter)
#Core_Mucus$aSlin <-predict(Adiameter)
### Anomalous Test
cor.test(log10(Core_Mucus$Diameter), Core_Mucus$alpha, 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test( log10(Core_Mucus$Diameter), Core_Mucus$alpha, 
          method = "pearson",use = "na.or.complete")
cor.test(Core_Mucus$Zeta,Core_Mucus$alpha, 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Core_Mucus$Zeta,Core_Mucus$alpha, 
         method = "pearson",use = "na.or.complete")

### Effective Diffusion
#Core_Mucus$Diffusion_constant = log10(Core_Mucus$Diffusion_constant)
#Core_Mucus$Diameter = log10(Core_Mucus$Diameter)

### Linear Analysis ###        
diameterlm = lm(log10(Diffusion_constant)~log10(Diameter), data = Core_Mucus)
Neglin = lm(log10(Diffusion_constant)~Zeta, data = Negative)
Poslin = lm(log10(Diffusion_constant)~Zeta, data = Positive)
Acharge = lm(alpha ~ Zeta, data = Core_Mucus)

summary(diameterlm)
summary(Neglin)
summary(Poslin)
summary(Acharge)

### Spearman and Pearson Analysis ###
###Diameter###
cor.test(log10(Core_Mucus$Diameter),log10(Core_Mucus$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Core_Mucus$Diameter),log10(Core_Mucus$Diffusion_constant), 
         method = "pearson",use = "na.or.complete")
### Negative Zeta ###
cor.test(Negative$Zeta,log10(Negative$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Negative$Zeta,log10(Negative$Diffusion_constant), 
         method = "pearson",use = "na.or.complete")
###Positive Zeta Potential###
cor.test(Positive$Zeta,log10(Positive$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Positive$Zeta,log10(Positive$Diffusion_constant), 
         method = "pearson",use = "na.or.complete")

# Old Plots ---------------------------------------------------------------


p1 <- ggplot(Core_Mucus, aes(x = Diameter, y = Diffusion_constant))
plt1 <- p1 + geom_point(size = .7, aes(shape = Surface_Chemistry ,color = alpha<0.8),stroke = .2)+
  scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) + scale_color_manual(na.value = "red", values = c("blue","black"), labels = c("~1.0","< 0.8")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  scale_y_log10(limits = c(1E-5,100), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10(limits = c(1,1500)) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"))+
  theme(axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4))) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Diameter, d (nm)")+theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank(), legend.background = element_rect(fill = "transparent"))+ theme(legend.text = element_text(size = 3.65))+ theme(legend.key.width = unit(.8,"line"))+
  theme(legend.spacing = unit(0,"cm")) + theme(legend.title = element_blank()) + theme(legend.justification = c(0,0),legend.position = c(.001,-.05))+ theme(panel.grid.major = element_blank())+
  geom_smooth(data = SmallMucus, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") + guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) +
  geom_vline(xintercept = 100, linetype = "dotted", size = 0.2)
#  
#
#
p2 <- ggplot(Core_Mucus, aes(x = Diameter, y = alpha))
plt2 <- p2 +  geom_point(size = .7, aes(shape = Surface_Chemistry),stroke = .2) + scale_shape_manual(values = c(6,4,2,1,0,5))+ theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4)) ) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha),x = "Diameter, d (nm)") + theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_log10(limits = c(1,1500))+ scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) +
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))


figlayout <- egg::ggarrange(plt1, plt2, nrow = 1)
#ggarrange(plots = list(p1,p2), nrow = 1, ncol = 2, widths = c(1,1))

ggsave("FIGUREP.pdf", figlayout, width = 10, height = 5, units = c("cm"), dpi = 2000, scale = 1:1)


c1 <- ggplot(Core_Mucus, aes(x = Zeta, y = Diffusion_constant))+geom_point(size = .7, aes(shape = Surface_Chemistry, color = alpha<0.8), stroke = .2)
c1plt <- c1 + scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses"))+ scale_color_manual(na.value = "red", values = c("blue","black"), labels = c("~1.0","< 0.8")) +
  scale_x_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80")) +
  scale_y_log10(limits = c(1E-5, 100), breaks = trans_breaks("log10",n= 8, function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+ theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4))) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Zeta-potential,"~zeta~"(mV)")+theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank())+ theme(legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"))+ theme(legend.key.width = unit(.3,"line"))+
  theme(legend.spacing = unit(0.4,"cm")) + theme(legend.title = element_blank()) + theme(legend.justification = c(0,0),legend.position = c(.02,.48))+ theme(panel.grid.major = element_blank())+
  geom_smooth(data = Negative, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") + 
  guides(color=guide_legend(override.aes=list(fill=NA)))+theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = "dotted", size =.2)
#geom_smooth(data = Positive, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "dashed")
c2 <- ggplot(Core_Mucus, aes(x = Zeta, y = alpha))+geom_point(size = .7, aes(shape = Surface_Chemistry),stroke = .2)+ scale_shape_manual(values = c(6,4,2,1,0,5))
c2plt <- c2 +
  theme(axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4)) ) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha), x = bquote("Zeta-potential,"~zeta~"(mV)")) + theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80,-60,-40,-20, 0, 20, 40, 60, 80) ) + scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) +theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(axis.ticks = element_line(size = rel(.3))) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

figlayout2 <- egg::ggarrange(c1plt, c2plt, nrow = 1)

ggsave("FIGUREC.pdf", figlayout2, width = 10, height = 5, units = c("cm"), dpi = 2000, scale = 1:1)



# Figure_2_Diffusion_versus_alpha -----------------------------------------



Diffvalpp <- ggplot(Core_Mucus, aes(x= alpha, y=Diffusion_constant))+
  geom_smooth(data = Core_Mucus, method = "lm", se = T, show.legend = F, level = 0.95, color = "black", size = .18, na.rm = TRUE)+ 
  geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2, color = "black", fill = "white")

Diffvalpplt <- Diffvalpp +  scale_shape_manual(values = c(25,4,24,21,22,23), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_log10(limits = c(1E-5, 1000),   breaks = c(1e-5,1e-3,1e-1,1e1,1e3), labels = trans_format("log10", math_format(10^.x))) +
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"), x = "Anomalous exponent,"~alpha) + 
  scale_x_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) +
  theme(axis.line = element_blank())+
  theme(legend.key = element_blank(), legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,0.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

t1 <- lm(log10(Diffusion_constant)~alpha, data = Core_Mucus)
summary(t1)

cor.test(Core_Mucus$alpha,log10(Core_Mucus$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Core_Mucus$alpha,log10(Core_Mucus$Diffusion_constant), 
         method = "pearson",use = "na.or.complete", exact = FALSE)
alpha_SLR <- lm(alpha~log10(Diffusion_constant), data = Core_Mucus)
summary(alpha_SLR)
predalb1 <- coefficients(alpha_SLR)
nonalpha_points <- data.frame(Diffusion_constant = NONALPHA$Diffusion_constant)
nonalpha_points_2 <- data.frame(Diffusion_constant = NONALPHA$Diffusion_constant, Surface_Chemistry = NONALPHA$Surface_Chemistry)
pred_alpha_points <- predict(alpha_SLR, nonalpha_points, interval = "prediction")
pred_alpha_points2 <- cbind(nonalpha_points_2, pred_alpha_points)
names(pred_alpha_points2)[3] <- "alpha"
pred_alpha_points2$alpha <- replace(pred_alpha_points2$alpha, pred_alpha_points2$alpha >= 1.0, 1.0 )
pred_alpha_points2$upr <- replace(pred_alpha_points2$upr, pred_alpha_points2$alpha >= 1.0, NA )
pred_alpha_points2$lwr <- replace(pred_alpha_points2$lwr, pred_alpha_points2$alpha >= 1.0, NA )

alphap <- ggplot(data = pred_alpha_points2, aes(x = Diffusion_constant, y = alpha)) + geom_smooth(data = Core_Mucus, method = "lm", se = T, show.legend = F, level = 0.95, color = "black", size = .18, na.rm = TRUE)+
  geom_point(size = .7, color = 'black', aes(shape = Surface_Chemistry), fill = "grey33",stroke = .2)
alphaplt <- alphap + scale_shape_manual(values = c(25,4,24,21,22,23), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) +
  
  scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) +
  scale_x_log10(limits = c(1E-5, 1000), breaks = c(1e-5,1e-3,1e-1,1e1,1e3), labels = trans_format("log10", math_format(10^.x))) + 
  theme(axis.text = element_text( size = rel(.4), colour = "black"), axis.ticks = element_line(size = rel(.3)),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) +
  labs(x = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"), y = "Anomalous exponent,"~alpha) +
  theme(legend.key = element_blank(), legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.45)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" )) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(axis.line = element_blank())+
  geom_line(aes(y = lwr), color = "black", linetype = 'dashed', size = .18) + geom_line(aes(y = upr), color = "black", linetype = 'dashed', size = .18) +
  geom_segment(aes(x = 3.6, xend = 100, y =1, yend = 1), size = .18, linetype = "solid") +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  annotation_logticks(side = "b", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 
#geom_vline(xintercept = 3.6, linetype = "dotted", size = 0.2)

DiffvAlpfiglayout <- egg::ggarrange(Diffvalpplt, alphaplt, nrow = 1)

ggsave("Figure_2_difvalp.pdf", DiffvAlpfiglayout, width = 12, height = 5.5, units = c("cm"), dpi = 2000)


# Labeling and Prediction -------------------------------------------------

Predexp <- data.frame(alpha = Core_Mucus$alpha)

Predexp[is.na(Predexp$alpha), "Data_Type"] <- "Prediction"

PREDT <- replace(Predexp$Data_Type, is.na(Predexp$Data_Type), "Experiment")
Dataclassification <- data.frame(Data_type = PREDT)

Mucus_Core1 <- replace(Core_Mucus$alpha, is.na(Core_Mucus$alpha), pred_alpha_points2$alpha)
Mucus_Core2 <- data.frame(alpha = Mucus_Core1)

removeMucus <- Core_Mucus
removeMucus$alpha <- NULL 
Mucus_Core3 <- cbind(removeMucus, Mucus_Core2)
Mucus_Core <- cbind(Mucus_Core3, Dataclassification)

SizeT1_ind     <- which(Mucus_Core$Diameter < 100.)
SizeT2_big_ind <- which(Mucus_Core$Diameter >= 100.)
alphasize_ind  <- which(Mucus_Core$alpha < 1)
alphasize1_ind <- which(Mucus_Core$alpha == 1)
positiveT_ind  <- which(Mucus_Core$Zeta > 0.0)
negativeT_ind  <- which(Mucus_Core$Zeta <= 0.0)

THCT_ind       <- which(Mucus_Core$Mucus_Source == "Human cervix" )
TPIT_ind       <- which(Mucus_Core$Mucus_Source == "Pig intestine")                                                 
TPST_ind       <- which(Mucus_Core$Mucus_Source == "Pig stomach")                                
THGT_ind      <- which(Mucus_Core$Mucus_Source == "Hydrogel")                                                                      
THLT_ind      <- which(Mucus_Core$Mucus_Source == "Human lung") 

Temp295_ind    <- which(Mucus_Core$Temperature == 295)
Temp298_ind    <- which(Mucus_Core$Temperature == 298)
Temp310_ind    <- which(Mucus_Core$Temperature == 310)  

pH3_ind       <-which((Mucus_Core$pH >= 3) &(Mucus_Core$pH < 4))
pH4_ind       <-which((Mucus_Core$pH >= 4) &(Mucus_Core$pH < 5))
pH6_ind       <-which((Mucus_Core$pH >= 6) &(Mucus_Core$pH < 7))
pH7_ind       <-which( (Mucus_Core$pH >= 7 ) & (Mucus_Core$pH <= 8) )

TMUC2_ind    <-which(Mucus_Core$Dominate_Mucin == "MUC2")
TMUC5AC_ind    <-which(Mucus_Core$Dominate_Mucin == "MUC5AC")
TMUC5B_ind    <-which(Mucus_Core$Dominate_Mucin == "MUC5B")

AlphaMucusT <- Mucus_Core[alphasize_ind, ]
AlphaMucus1T <- Mucus_Core[alphasize1_ind, ]
BigMucusT <- Mucus_Core[SizeT2_big_ind, ]
SmallMucusT <- Mucus_Core[SizeT1_ind, ]
PositiveT <- Mucus_Core[positiveT_ind, ]
NegativeT <- Mucus_Core[negativeT_ind, ]

TCervix   <- Mucus_Core[THCT_ind, ]
TIntestine <- Mucus_Core[TPIT_ind, ]
TStomach    <- Mucus_Core[TPST_ind, ]
THydrogel      <- Mucus_Core[THGT_ind, ]
TLung        <- Mucus_Core[THLT_ind, ]

Temp295T       <- Mucus_Core[Temp295_ind, ]
Temp298T       <- Mucus_Core[Temp298_ind, ]
Temp310T       <- Mucus_Core[Temp310_ind, ]

pH3T           <- Mucus_Core[pH3_ind, ]
pH4T           <- Mucus_Core[pH4_ind, ]
pH6T           <- Mucus_Core[pH6_ind, ]
pH7T           <- Mucus_Core[pH7_ind, ]

TMUC2    <- Mucus_Core[TMUC2_ind, ] 
TMUC5AC  <- Mucus_Core[TMUC5AC_ind, ]
TMUC5B   <- Mucus_Core[TMUC5B_ind, ]

# Mechanistic model lm ----------------------------------------------------
mechvec1 <- data.frame(Diffusion_constant = Core_Mucus$Diffusion_constant, alpha = 1-Core_Mucus$alpha, D_w = Core_Mucus$D_w)
mechvec <- na.omit(mechvec1)
MechModlm <- lm(log10(Diffusion_constant)~alpha, data = mechvec)
summary(MechModlm)
mean(mechvec$D_w)
sd(mechvec$D_w)


# Table_S3_linear_analysis ------------------------------------------------

sink("Table_S3_linear_analysis.txt")

# Linear analysis of effective diffusion versus diameter smaller than 100 nm
TSmallmucuslm = lm(log10(Diffusion_constant)~log10(Diameter), data = SmallMucusT)
TSmallmucus_spearman = cor.test(log10(SmallMucusT$Diameter),log10(SmallMucusT$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
TSmallmucus_pearson = cor.test(log10(SmallMucusT$Diameter),log10(SmallMucusT$Diffusion_constant),  
         method = "pearson",use = "na.or.complete")
summary(TSmallmucuslm)
TSmallmucus_spearman
TSmallmucus_pearson

# Linear analysis of effective diffusion versus negative zeta potential
TNeglin       = lm(log10(Diffusion_constant)~Zeta, data = NegativeT)
TNeglin_spearman = cor.test(NegativeT$Zeta,log10(NegativeT$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
TNeglin_pearson = cor.test(NegativeT$Zeta,log10(NegativeT$Diffusion_constant), 
         method = "pearson",use = "na.or.complete")
summary(TNeglin)
TNeglin_spearman
TNeglin_pearson

# Linear analysis of effective diffusion versus positive zeta potential
TPoslin       = lm(log10(Diffusion_constant)~Zeta, data = PositiveT)
TPoslin_spearman = cor.test(PositiveT$Zeta,log10(PositiveT$Diffusion_constant),
         method = "spearman",use = "na.or.complete", exact = FALSE)
TPoslin_pearson = cor.test(PositiveT$Zeta,log10(PositiveT$Diffusion_constant), 
         method = "pearson",use = "na.or.complete")
summary(TPoslin)
TPoslin_spearman
TPoslin_pearson

# Linear analysis of effective diffusion versus anomalous exponent
Figure_2_difvalp_lm <- lm(log10(Diffusion_constant)~alpha, data = Core_Mucus)
Figure_2_difvalp_spearman = cor.test(Core_Mucus$alpha,log10(Core_Mucus$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
Figure_2_difvalp_pearson = cor.test(Core_Mucus$alpha,log10(Core_Mucus$Diffusion_constant),  
         method = "pearson",use = "na.or.complete")
summary(Figure_2_difvalp_lm)
Figure_2_difvalp_spearman
Figure_2_difvalp_pearson

# Linear analysis of anomalous exponent versus diameter less than 100 nm
AsmalldiaT1 = lm(alpha ~ log10(Diameter), data = SmallMucusT)
AsmalldiaT1_spearman = cor.test(log10(SmallMucusT$Diameter),SmallMucusT$alpha, 
                                method = "spearman",use = "na.or.complete", exact = FALSE)
AsmalldiaT1_pearson = cor.test(log10(SmallMucusT$Diameter),SmallMucusT$alpha, 
                               method = "pearson",use = "na.or.complete")
summary(AsmalldiaT1)
AsmalldiaT1_spearman
AsmalldiaT1_pearson

# Linear analysis of anomalous exponent versus netaive zeta potential 
TANeglin      = lm(alpha ~ Zeta, data = NegativeT)
TANeglin_spearman = cor.test(NegativeT$Zeta, NegativeT$alpha, 
                             method = "spearman",use = "na.or.complete", exact = FALSE)
TANeglin_pearson = cor.test(NegativeT$Zeta, NegativeT$alpha, 
                            method = "pearson",use = "na.or.complete")
summary(TANeglin)
TANeglin_spearman
TANeglin_pearson

# Linear analysis of anomalous exponentversus positive zeta potential 
TAPoslin      = lm(alpha ~ Zeta, data = PositiveT)
TAPoslin_spearman = cor.test(PositiveT$Zeta, PositiveT$alpha, 
                             method = "spearman",use = "na.or.complete", exact = FALSE)
TAPoslin_pearson = cor.test(PositiveT$Zeta, PositiveT$alpha, 
                            method = "pearson",use = "na.or.complete")
summary(TAPoslin)
TAPoslin_spearman
TAPoslin_pearson

sink()
# Data files for linear analysis 
Table_S3_a <- data.frame(Diameter = SmallMucusT$Diameter, Diffusion_constant = SmallMucusT$Diffusion_constant, Particle_type = SmallMucusT$Surface_Chemistry)
Table_S3_b <- data.frame(Zeta = NegativeT$Zeta, Diffusion_constant = NegativeT$Diffusion_constant            , Particle_type = NegativeT$Surface_Chemistry)
Table_S3_c <- data.frame(Zeta = PositiveT$Zeta, Diffusion_constant = PositiveT$Diffusion_constant            , Particle_type = PositiveT$Surface_Chemistry)
Table_S3_d <- data.frame(alpha = Core_Mucus$alpha, Diffusion_constant = Core_Mucus$Diffusion_constant        , Particle_type = Core_Mucus$Surface_Chemistry)
Table_S3_e <- data.frame(Diameter = AlphaMucusT$Diameter, alpha = AlphaMucusT$alpha                          , Particle_type = AlphaMucusT$Surface_Chemistry)
Table_S3_f <- data.frame(Zeta = NegativeT$Zeta, alpha = NegativeT$alpha                                      , Particle_type = NegativeT$Surface_Chemistry)
Table_S3_g <- data.frame(Zeta = PositiveT$Zeta, alpha = PositiveT$alpha                                      , Particle_type = PositiveT$Surface_Chemistry)
    
write.csv(Table_S3_a, "Table_S3_a.csv", row.names = FALSE)
write.csv(Table_S3_b, "Table_S3_b.csv", row.names = FALSE)
write.csv(Table_S3_c, "Table_S3_c.csv", row.names = FALSE)
write.csv(Table_S3_d, "Table_S3_d.csv", row.names = FALSE)
write.csv(Table_S3_e, "Table_S3_e.csv", row.names = FALSE)
write.csv(Table_S3_f, "Table_S3_f.csv", row.names = FALSE)
write.csv(Table_S3_g, "Table_S3_g.csv", row.names = FALSE)



# Additional_linear_analysis ----------------------------------------------


AdiameterT1 = lm(alpha ~ log10(Diameter), data = AlphaMucusT)



AdiameterT2 = lm(alpha ~ log10(Diameter), data = Mucus_Core)
summary(AdiameterT2)


ATemperatureT = lm(alpha ~ log10(Temperature), data = Mucus_Core)
cor.test(log10(Mucus_Core$Temperature),Mucus_Core$alpha, 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Mucus_Core$Temperature),Mucus_Core$alpha, 
         method = "pearson",use = "na.or.complete")
summary(ATemperatureT)


ApHT = lm(alpha ~ pH, data = Mucus_Core)
cor.test(Mucus_Core$pH,Mucus_Core$alpha, 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Mucus_Core$pH,Mucus_Core$alpha, 
         method = "pearson",use = "na.or.complete")
summary(ApHT)





Tdiameterlm   = lm(log10(Diffusion_constant)~log10(Diameter), data = Mucus_Core)
cor.test(log10(Mucus_Core$Diameter),log10(Mucus_Core$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Mucus_Core$Diameter),log10(Mucus_Core$Diffusion_constant),  
         method = "pearson",use = "na.or.complete")
summary(Tdiameterlm)


Tdiasmalllm   = lm(log10(Diffusion_constant)~log10(Diameter), data = AlphaMucus1T)
summary(Tdiasmalllm)

Tcharlin      = lm(log10(Diffusion_constant)~Zeta, data = Mucus_Core)
cor.test(Mucus_Core$Zeta,log10(Mucus_Core$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Mucus_Core$Zeta,log10(Mucus_Core$Diffusion_constant), 
         method = "pearson",use = "na.or.complete")
summary(Tcharlin)


TAcharge      = lm(alpha ~ Zeta, data = Mucus_Core)
summary(TAcharge)


TTemperaturelm = lm(log10(Diffusion_constant)~log10(Temperature), data = Mucus_Core)
cor.test(log10(Mucus_Core$Temperature),log10(Mucus_Core$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Mucus_Core$Temperature),log10(Mucus_Core$Diffusion_constant),  
         method = "pearson",use = "na.or.complete")
summary(TTemperaturelm)


TpHlm = lm(log10(Diffusion_constant)~pH, data = Mucus_Core)
cor.test(Mucus_Core$pH,log10(Mucus_Core$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Mucus_Core$pH,log10(Mucus_Core$Diffusion_constant),  
         method = "pearson",use = "na.or.complete")

summary(TpHlm)




# Figure_4_Particle_size --------------------------------------------------


Tp1 <- ggplot(Mucus_Core, aes(x = Diameter, y = Diffusion_constant)) + 
  geom_smooth(data = AlphaMucus1T, method = "lm", se = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "dashed") +
  geom_smooth(data = SmallMucusT, method = "lm", se = T,fullrange = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid")
Tplt1 <- Tp1 + geom_point(size = .7, aes(shape = Surface_Chemistry,color = Data_type, fill = Data_type),stroke = .2)+
  scale_shape_manual(values = c(25,4,22,21,24,23), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) +
  scale_color_manual(values = c("Black", "Black"), labels = c("Empirical data", "Predicted data")) +
  scale_fill_manual(values = c("white", "white"), labels = c("Empirical data", "Predicted data")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  scale_y_log10(limits = c(1E-5,1000), breaks = c(1e-5,1e-3,1e-1,1e1,1e3),  labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Diameter, d (nm)")+
  theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank(), legend.background = element_rect(fill = "transparent"))+ 
  theme(legend.text = element_text(size = rel(.35)))+ 
  theme(legend.key.width = unit(.8,"line"))+
  theme(legend.spacing = unit(0,"cm")) + 
  theme(legend.title = element_blank()) + 
  theme(legend.justification = c(0,0),legend.position = c(.001,-.05))+ 
  theme(panel.grid.major = element_blank())+
  theme(axis.line = element_blank())+
  guides(color=guide_legend(override.aes=list(fill=NA)))+ 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  theme(panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 100, linetype = "dotted", size = 0.2) + 
  annotation_logticks(side = "lb", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


Tp2 <- ggplot(Mucus_Core, aes(x = Diameter, y = alpha))+ geom_smooth(data = SmallMucusT, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid")
Tplt2 <- Tp2 +  geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type, fill = Data_type), stroke = .2) +
  scale_shape_manual(values = c(25,4,22,21,24,23))+ 
  scale_color_manual(values = c("Black", "Black"), labels = c("Empirical data", "Predicted data")) +
  scale_fill_manual(values = c("white", "grey33"), labels = c("Empirical data", "Predicted data")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha),x = "Diameter, d (nm)") + 
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(limits = c(0,1.5), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) + 
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.line = element_blank())+
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  coord_cartesian(xlim = c(1,1500), ylim = c(0,1.1), clip = "off")+
  geom_vline(xintercept = 100, linetype = "dotted", size = 0.2) + border("black", size = .65) +
  geom_segment(aes(x=3.5, xend = 55, y=1.0, yend = 1.0), size = .18, linetype = "dashed") + 
  annotation_logticks(side = "b", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) 



Tfiglayout <- egg::ggarrange(Tplt1, Tplt2, nrow = 1)


ggsave("Figure_4_Size.pdf", Tfiglayout, width = 12, height = 5.5, units = c("cm"), dpi = 2000, scale = 1:1)




# Reviewer's comments -----------------------------------------------------



Tp1s <- ggplot(Mucus_Core, aes(x = Diameter*1e-3, y = Diffusion_constant*Diameter*1e-3)) + 
  geom_smooth(data = AlphaMucus1T, method = "lm", se = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "dashed") +
  geom_smooth(data = SmallMucusT, method = "lm", se = T,fullrange = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid")
Tplt1s <- Tp1s + geom_point(size = .7, aes(shape = Surface_Chemistry,color = Data_type, fill = Data_type),stroke = .2)+
  scale_shape_manual(values = c(25,4,22,21,24,23), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) +
  scale_color_manual(values = c("Black", "Black"), labels = c("Empirical data", "Predicted data")) +
  scale_fill_manual(values = c("white", "white"), labels = c("Empirical data", "Predicted data")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  scale_y_log10(limits = c(1e-6,1e1), breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10(limits = c(1e-3,1e1),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote(""~D[eff]~"*d ("~mu~m^3~"/s)"),x = "d ("~mu~"m)")+
  theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank(), legend.background = element_rect(fill = "transparent"))+ 
  theme(legend.text = element_text(size = rel(.35)))+ 
  theme(legend.key.width = unit(.8,"line"))+
  theme(legend.spacing = unit(0,"cm")) + 
  theme(legend.title = element_blank()) + 
  theme(legend.justification = c(0,0),legend.position = c(.001,-.05))+ 
  theme(panel.grid.major = element_blank())+
  theme(axis.line = element_blank())+
  guides(color=guide_legend(override.aes=list(fill=NA)))+ 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  theme(panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 100*1e-3, linetype = "dotted", size = 0.2) + 
  annotation_logticks(side = "lb", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


Tp2s <- ggplot(Mucus_Core, aes(x = Diameter*1e-3, y = alpha*Diameter*1e-3))+ 
  geom_smooth(data = SmallMucusT, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid")+
  geom_smooth(data = AlphaMucus1T, method = "lm", se = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "dashed")+
  geom_smooth(data = BigMucusT, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "twodash")
Tplt2s <- Tp2s +  geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type, fill = Data_type), stroke = .2) +
  scale_shape_manual(values = c(25,4,22,21,24,23))+ 
  scale_color_manual(values = c("Black", "Black"), labels = c("Empirical data", "Predicted data")) +
  scale_fill_manual(values = c("white", "grey33"), labels = c("Empirical data", "Predicted data")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote(""~alpha~"*d ("~mu~"m)"),x = "d ("~mu~"m)") + 
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_log10(limits = c(1e-3,1e1),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits = c(1e-6, 1e1),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) + 
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.line = element_blank())+
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  coord_cartesian(clip = "off")+
  geom_vline(xintercept = 100*1e-3, linetype = "dotted", size = 0.2) + border("black", size = .65) +
  annotation_logticks(side = "lb", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) 



Tfiglayouts2 <- egg::ggarrange(Tplt1s, Tplt2s, nrow = 1)


ggsave("Figure_S1_diffuion_scaled_by_size.pdf", Tfiglayouts2, width = 12, height = 5.5, units = c("cm"), dpi = 2000, scale = 1:1)


ggsave("Figure_S1_diffuion_scaled_by_size_a.pdf", Tplt1s, width = 12, height = 5.5, units = c("cm"), dpi = 2000, scale = 1:1)

Table_Fig_S1_a <- data.frame(Zeta = Mucus_Core$Diameter, Diffusion_constant = Mucus_Core$Diffusion_constant                                      , Particle_type = Mucus_Core$Surface_Chemistry)

write.csv(Table_Fig_S1_a, "Figure_S1_diffuion_scaled_by_size_data_a.csv", row.names = FALSE)


#sink("ResearchCodeReviewer.txt")


TdiameterlmDD   = lm(log10(Diffusion_constant*Diameter*1e-3)~log10(Diameter*1e-3), data = Mucus_Core)
cor.test(log10(Mucus_Core$Diameter*1e-3),log10(Mucus_Core$Diffusion_constant*Mucus_Core$Diameter*1e-3), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Mucus_Core$Diameter*1e-3),log10(Mucus_Core$Diffusion_constant*Mucus_Core$Diameter*1e-3),  
         method = "pearson",use = "na.or.complete")
summary(TdiameterlmDD)


TdiameterlmDDsmall   = lm(log10(Diffusion_constant*Diameter*1e-3)~log10(Diameter*1e-3), data = SmallMucusT)
cor.test(log10(Mucus_Core$Diameter*1e-3),log10(Mucus_Core$Diffusion_constant*Mucus_Core$Diameter*1e-3), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Mucus_Core$Diameter*1e-3),log10(Mucus_Core$Diffusion_constant*Mucus_Core$Diameter*1e-3),  
         method = "pearson",use = "na.or.complete")
summary(TdiameterlmDDsmall)



TdiameterlmAD   = lm(log10(Diffusion_constant*Diameter*1e-3)~log10(Diameter*1e-3), data = Mucus_Core)
cor.test(log10(Mucus_Core$Diameter*1e-3),log10(Mucus_Core$alpha*Mucus_Core$Diameter*1e-3), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(Mucus_Core$Diameter*1e-3),log10(Mucus_Core$alpha*Mucus_Core$Diameter*1e-3),  
         method = "pearson",use = "na.or.complete")
summary(TdiameterlmAD)


TdiameterlmADsmall   = lm(log10(alpha*Diameter*1e-3)~log10(Diameter*1e-3), data = SmallMucusT)
cor.test(log10(SmallMucusT$Diameter*1e-3),log10(SmallMucusT$alpha*SmallMucusT$Diameter*1e-3), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(SmallMucusT$Diameter*1e-3),log10(SmallMucusT$alpha*SmallMucusT$Diameter*1e-3),  
         method = "pearson",use = "na.or.complete")
summary(TdiameterlmADsmall)



TdiameterlmADbig   = lm(log10(alpha*Diameter*1e-3)~log10(Diameter*1e-3), data = BigMucusT)
cor.test(log10(BigMucusT$Diameter*1e-3),log10(BigMucusT$alpha*BigMucusT$Diameter*1e-3), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(log10(BigMucusT$Diameter*1e-3),log10(BigMucusT$alpha*BigMucusT$Diameter*1e-3),  
         method = "pearson",use = "na.or.complete")
summary(TdiameterlmADbig)

#sink()

# Figure_3_Net_Charge -----------------------------------------------------



Tc1 <- ggplot(Mucus_Core, aes(x = Zeta, y = Diffusion_constant)) + geom_smooth(data = NegativeT, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") 

Tc1plt <- Tc1 + geom_point(size = .7, aes(shape = Surface_Chemistry,color = Data_type, fill = Data_type), stroke = .2) +
  scale_shape_manual(values = c(25,4,22,21,24,23), labels = c("AMINE","Antibodies              
 and proteins","Chitosan","COOH", "PEG","Viruses"))+
  scale_color_manual(values = c("Black", "Black"), labels = c("Empirical data", "Predicted data")) +
  scale_fill_manual(values = c("white", "white"), labels = c("Empirical data", "Predicted data")) +
  scale_x_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80")) +
  scale_y_log10(limits = c(1E-5, 1000), breaks = c(1e-5,1e-3,1e-1,1e1,1e3),  labels = trans_format("log10", math_format(10^.x)))  +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+ 
  theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) + 
  theme(axis.line = element_blank())+
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Zeta-potential,"~zeta~"(mV)")+
  theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank())+ 
  theme(legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"))+ 
  theme(legend.key.width = unit(.3,"line"))+
  theme(legend.spacing = unit(0.4,"cm")) + 
  theme(legend.title = element_blank()) +
  theme(legend.justification = c(0,0),legend.position = c(.02,.41))+ 
  theme(panel.grid.major = element_blank())+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.grid.minor = element_blank())+ 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  geom_vline(xintercept = 0, linetype = "dotted", size =.2) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


Tc2 <- ggplot(Mucus_Core, aes(x = Zeta, y = alpha))+  
  geom_smooth(data = NegativeT, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid")+
  geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type, fill = Data_type), stroke = .2)+
  scale_shape_manual(values = c(25,4,22,21,24,23))
Tc2plt <- Tc2 + scale_color_manual(values = c("Black", "black"), labels = c("Empirical data", "Predicted data")) +
  scale_fill_manual(values = c("white", "grey33"), labels = c("Empirical data", "Predicted data")) +
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) +
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha), x = bquote("Zeta-potential,"~zeta~"(mV)")) + 
  theme(legend.position = "none") + theme(axis.line = element_blank())+
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80,-60,-40,-20, 0, 20, 40, 60, 80) ) + 
  geom_vline(xintercept = 0, linetype = "dotted", size =.2) + 
  scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) +
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme(axis.ticks = element_line(size = rel(.3))) + 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

Tfiglayout2 <- egg::ggarrange(Tc1plt, Tc2plt, nrow = 1)

ggsave("Figure_3_Charge.pdf", Tfiglayout2, width = 12, height = 5.5, units = c("cm"), dpi = 2000, scale = 1:1)


# Temperature -------------------------------------------------------------



Temp295onD <- data.frame(Temp_group = "295 K",  Diffusion_constant = Temp295T$Diffusion_constant)
Temp298onD <- data.frame(Temp_group = "298 K",  Diffusion_constant = Temp298T$Diffusion_constant)
Temp310onD <- data.frame(Temp_group = "310 K",  Diffusion_constant = Temp310T$Diffusion_constant)
Temp295onA <- data.frame(Temp_group = "295 K",  alpha = Temp295T$alpha)
Temp298onA <- data.frame(Temp_group = "298 K",  alpha = Temp298T$alpha)
Temp310onA <- data.frame(Temp_group = "310 K",  alpha = Temp310T$alpha)

TempTboxD <- rbind(Temp295onD,Temp298onD,Temp310onD)
TempTboxA <- rbind(Temp295onA,Temp298onA,Temp310onA)


# pH ----------------------------------------------------------------------



TpH1 <- ggplot(Mucus_Core, aes(x = pH, y = Diffusion_constant))

TpH1plt <- TpH1 + geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2) +
  scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies              
                                                         and Proteins","Chitosan","COOH", "PEG","Viruses"))+
  scale_x_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0")) +
  scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))  +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+ 
  theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) +
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "pH")+theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank())+ 
  theme(legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"))+
  theme(legend.key.width = unit(.3,"line"))+
  theme(legend.spacing = unit(0.4,"cm")) + 
  theme(legend.title = element_blank()) + 
  theme(legend.justification = c(0,0),legend.position = c(.02,.48))+
  theme(panel.grid.major = element_blank())+
  guides(color=guide_legend(override.aes=list(fill=NA)))+ 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  theme(panel.grid.minor = element_blank())

TpH2 <- ggplot(Mucus_Core, aes(x = pH, y = alpha))+geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type),stroke = .2)+ scale_shape_manual(values = c(6,4,2,1,0,5))
TpH2plt <- TpH2 + scale_color_manual(values = c("Black", "lightsteelblue3"), labels = c("Experiment", "Prediction")) +
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha), x = bquote("pH")) + theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0")) +  
  scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) +theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(axis.ticks = element_line(size = rel(.3))) +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

pH3onD <- data.frame(pH_group = "pH 3",  Diffusion_constant = pH3T$Diffusion_constant)
pH4onD <- data.frame(pH_group = "pH 4",  Diffusion_constant = pH4T$Diffusion_constant)
pH6onD <- data.frame(pH_group = "pH ~6.5",  Diffusion_constant = pH6T$Diffusion_constant)
pH7onD <- data.frame(pH_group = "pH ~7", Diffusion_constant = pH7T$Diffusion_constant)
pH3onA <- data.frame(pH_group = "pH 3",  alpha = pH3T$alpha)
pH4onA <- data.frame(pH_group = "pH 4",  alpha = pH4T$alpha)
pH6onA <- data.frame(pH_group = "pH ~6.5",  alpha = pH6T$alpha)
pH7onA <- data.frame(pH_group = "pH ~7", alpha = pH7T$alpha)

pHTotD <- rbind(pH3onD,pH4onD,pH6onD,pH7onD)
pHTotA <- rbind(pH3onA,pH4onA,pH6onA,pH7onA)

TpHboxDIF <- ggboxplot(pHTotD, x = "pH_group", y = "Diffusion_constant",width = .5, lwd = .2)
TpHboxpltDIF  <- TpHboxDIF+ scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Effective diffusion,"~D[eff]~"("~mu~m^2~"/s)"))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
TpHboxALP <- ggboxplot(pHTotA, x = "pH_group", y = "alpha",width = .5, lwd = .2)
TpHboxpltALP <- TpHboxALP + scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.50,.75,1.0), labels = c("0.00", "0.25", "0.50", "0.75","1.00")) +
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text(size = rel(.4), color = "Black"), 
        axis.title.x = element_text(size = rel(.45)), 
        axis.title.y = element_text(size = rel(.45)))+
  theme(axis.ticks = element_line(size = rel(.3)))

TfiglayoutpH <- egg::ggarrange(TpH1plt, TpH2plt, TpHboxpltDIF, TpHboxpltALP, nrow = 2, ncol = 2)

ggsave("FIGUREpH.pdf", TfiglayoutpH, width = 14, height = 12, units = c("cm"), dpi = 2000, scale = 1:1)


# Mucus and Salt Con ------------------------------------------------------



MUCCONdiflm <- lm(log10(Diffusion_constant)~Muc_Con, data = Mucus_Core)
summary(MUCCONdiflm)
MUCCONalplm <- lm(alpha~Muc_Con, data = Mucus_Core)
summary(MUCCONalplm)
SALCONdiflm <- lm(log10(Diffusion_constant)~Salt_Concentration, data = Mucus_Core)
summary(SALCONdiflm)
SALCONalplm <- lm(alpha~Salt_Concentration, data = Mucus_Core)
summary(SALCONalplm)


TMCC1 <- ggplot(Mucus_Core, aes(x = Muc_Con, y = Diffusion_constant)) 
TMUCCON <- TMCC1 + geom_point(size = .7, aes(shape = Surface_Chemistry),stroke = .2)+
  scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
                                                         and proteins","Chitosan","COOH", "PEG","Viruses")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  scale_y_log10(limits = c(1E-5,1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(0,5),breaks = c(0,1,2,3,4),  labels = c("0","1","2","3","4")) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Mucus concentration (% wt/vol)")+theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank(), legend.background = element_rect(fill = "transparent"))+
  theme(legend.text = element_text(size = rel(.45)))+ theme(legend.key.width = unit(.8,"line"))+
  theme(legend.spacing = unit(0,"cm")) + theme(legend.title = element_blank()) + theme(legend.justification = c(0,0),legend.position = c(.51,0))+ theme(panel.grid.major = element_blank())+
  geom_smooth(data = Mucus_Core, method = "lm", se = T,fullrange = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") +
  guides(color=guide_legend(override.aes=list(fill=NA)))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  theme(panel.grid.minor = element_blank())

TMCC2 <- ggplot(Mucus_Core, aes(x = Muc_Con, y = alpha)) 
TMCCplt2 <- TMCC2 +  geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type),stroke = .2) + scale_shape_manual(values = c(6,4,2,1,0,5))+ theme(axis.ticks = element_line(size = rel(.3)))+
  scale_color_manual(values = c("Black", "lightsteelblue3"), labels = c("Experiment", "Prediction")) +
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha),x = "Mucus concentration (% wt/vol)") + 
  theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(0,5),breaks = c(0,1,2,3,4),  labels = c("0","1","2","3","4"))+
  scale_y_continuous(limits = c(0,1.5), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) + theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) + coord_cartesian(xlim = c(0,5), ylim = c(0,1.1))+
  geom_smooth(data = Mucus_Core, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) 

TSCC1 <- ggplot(Mucus_Core, aes(x = Salt_Concentration, y = Diffusion_constant)) 
TSALCON <- TSCC1 + geom_point(size = .7, aes(shape = Surface_Chemistry),stroke = .2)+
  scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
                                                         and proteins","Chitosan","COOH", "PEG","Viruses")) +
  theme(axis.ticks = element_line(size = rel(.3)))+
  scale_y_log10(limits = c(1E-5,1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(0,510),breaks = c(0,100,200,300,400,500),  labels = c("0","100","200","300","400","500")) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Salt concentration (mM)")+theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank(), legend.background = element_rect(fill = "transparent"))+ 
  theme(legend.text = element_text(size = rel(.45)))+ theme(legend.key.width = unit(.8,"line"))+
  theme(legend.spacing = unit(0,"cm")) + theme(legend.title = element_blank()) + theme(legend.justification = c(0,0),legend.position = c(.51,0))+ theme(panel.grid.major = element_blank())+
  geom_smooth(data = Mucus_Core, method = "lm", se = T,fullrange = F, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") +
  guides(color=guide_legend(override.aes=list(fill=NA)))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  theme(panel.grid.minor = element_blank())

TSCC2 <- ggplot(Mucus_Core, aes(x = Salt_Concentration, y = alpha)) 
TSCCplt2 <- TSCC2 +  geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type),stroke = .2) + scale_shape_manual(values = c(6,4,2,1,0,5))+ theme(axis.ticks = element_line(size = rel(.3)))+
  scale_color_manual(values = c("Black", "lightsteelblue3"), labels = c("Experiment", "Prediction")) +
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha),x = "Salt concentration (mM)") + theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(0,510),breaks = c(0,100,200,300,400,500),  labels = c("0","100","200","300","400","500"))+
  scale_y_continuous(limits = c(0,1.5), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) + theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) + coord_cartesian(xlim = c(0,510), ylim = c(0,1.1))+
  geom_smooth(data = Mucus_Core, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

Concenfigurelayout<-egg::ggarrange(TMUCCON,TMCCplt2,TSALCON,TSCCplt2, nrow = 2, ncol = 2)
ggsave("FIGUREconcentration.pdf", Concenfigurelayout, width = 14, height = 12, units = c("cm"), dpi = 2000, scale = 1:1)


# Figure_S4_Indepth_SLR_analysis ------------------------------------------





TRDA   <- lm(log10(Diffusion_constant)~alpha, data = Core_Mucus, na.action = na.exclude)
TRDDS  <- lm(log10(Diffusion_constant)~log10(Diameter), data = SmallMucusT, na.action = na.exclude)   # log10(Deff) vs log10(diameter) for d =< 100
TRDDA  <- lm(log10(Diffusion_constant)~log10(Diameter), data = AlphaMucus1T, na.action = na.exclude)  # log10(Deff) vs Diameter for alpha = 1
TRDZ   <- lm(log10(Diffusion_constant)~Zeta, data = NegativeT, na.action = na.exclude) # log10(Deff) vs negative zeta

RES1 <- ggplot(fortify(TRDA, Core_Mucus), aes(x= alpha, .resid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
RES1plt <- RES1 + scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
                                                                     and proteins","Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion residuals"), x = "Anomalous exponent,"~alpha) + 
  scale_x_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) +
  theme(axis.line = element_blank()) +
  theme(axis.text = element_text(size = rel(1), colour = "black"),
        axis.ticks = element_line(size = rel(.3)),
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.title.x = element_text(size = rel(.45)), 
        axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), 
        legend.text = element_text(size = rel(.35)),
        legend.background = element_rect(fill = "transparent"), 
        legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),
        legend.title = element_blank(),
        legend.justification = c(0,0),
        legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white" ,colour = "black"))+ 
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2)+
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) 

RES2 <- ggplot(fortify(TRDDS, SmallMucusT), aes(x= Diameter, .resid))+
  geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
RES2plt <- RES2 + scale_shape_manual(values = c(4,2,1,0,5), labels = c("Antibodies 
                                                                     and proteins","Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion residuals"), x = "Diameter, d (nm)") + 
  scale_x_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text( size = rel(1), colour = "black"),
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)),
        axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(),
        legend.text = element_text(size = rel(.35)),
        legend.background = element_rect(fill = "transparent"),
        legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), 
        legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white" ,colour = "black")) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed",size = .2) + 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  geom_vline(xintercept = 100, linetype = "dotted", size =.2) + 
  annotation_logticks(side = "b", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


RES3 <- ggplot(fortify(TRDDA, AlphaMucus1T), aes(x= Diameter, .resid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
RES3plt <- RES3 + scale_shape_manual(values = c(4,5), labels = c("Antibodies 
and proteins","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion residuals"), x = "Diameter, d (nm)") + 
  scale_x_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), axis.line = element_blank(),
        axis.ticks = element_line(size = rel(.3)),
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.title.x = element_text(size = rel(.45)),
        axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) +
  geom_vline(xintercept = 100, linetype = "dotted", size =.2) + 
  annotation_logticks(side = "b", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


RES4 <- ggplot(fortify(TRDZ, NegativeT), aes(x= Zeta, y= .resid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
RES4plt <- RES4 + scale_shape_manual(values = c(2,1,0,5), labels = c("Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion residuals"), x = "Zeta-potential,"~zeta~"(mV)") + 
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80,-60,-40,-20, 0, 20, 40, 60, 80) ) + 
  geom_vline(xintercept = 0, linetype = "dotted", size =.2) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

#### Standardized Residuals
SRES1 <- ggplot(fortify(TRDA, Core_Mucus), aes(x= alpha, .stdresid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
SRES1plt <- SRES1 + scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
                                                                         and proteins","Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Anomalous exponent,"~alpha) + 
  scale_x_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00")) +
  theme(axis.text = element_text( size = rel(1), colour = "black"),
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = rel(.35)),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

SRES2 <- ggplot(fortify(TRDDS, SmallMucusT), aes(x= Diameter, .stdresid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
SRES2plt <- SRES2 + scale_shape_manual(values = c(4,2,1,0,5), labels = c("Antibodies 
                                                                     and proteins","Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Diameter, d (nm)") + 
  scale_x_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2) +theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  geom_vline(xintercept = 100, linetype = "dotted", size =.2) + 
  annotation_logticks(side = "b", size = 0.185, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


SRES3 <- ggplot(fortify(TRDDA, AlphaMucus1T), aes(x= Diameter, .stdresid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
SRES3plt <- SRES3 + scale_shape_manual(values = c(4,5), labels = c("Antibodies 
and proteins","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Diameter, d (nm)") + 
  scale_x_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line= element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  geom_vline(xintercept = 100, linetype = "dotted", size =.2) + 
  annotation_logticks(side = "b", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 


SRES4 <- ggplot(fortify(TRDZ, NegativeT), aes(x= Zeta, .stdresid))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
SRES4plt <- SRES4 + scale_shape_manual(values = c(2,1,0,5), labels = c("Chitosan","COOH", "PEG","Viruses")) + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = bquote("Zeta-potential,"~zeta~"(mV)")) + 
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80,-60,-40,-20, 0, 20, 40, 60, 80) ) + 
  geom_vline(xintercept = 0, linetype = "dotted", size =.2) +
  theme(axis.text = element_text( size = rel(1), colour = "black"),
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .2)+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

#### Density plots for Residuals

DRES1 <- ggplot(fortify(TRDA, Core_Mucus))+ stat_qq(size = .7, aes(sample = .stdresid), stroke = .2) + stat_qq_line(aes(sample = .stdresid), size = 0.18) 
DRES1plt <- DRES1 + scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Theoretical quantiles of Anomalous exponent") + 
  scale_x_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3","-2", "-1", "0", "1", "2","3")) +
  theme(axis.text = element_text( size = rel(1), colour = "black"),
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

DRES2 <- ggplot(fortify(TRDDS, SmallMucusT))+ stat_qq(size = .7, aes(sample = .stdresid), stroke = .2) + stat_qq_line(aes(sample = .stdresid),size = 0.18)
DRES2plt <- DRES2 + scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Theorectical quantiles of diameter") + 
  scale_x_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3","-2", "-1", "0", "1", "2","3")) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + 
  guides(color=guide_legend(override.aes=list(fill=NA)))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))

DRES3 <- ggplot(fortify(TRDDA, AlphaMucus1T))+ stat_qq(size = .7, aes(sample = .stdresid), stroke = .2) + stat_qq_line(aes(sample = .stdresid),size = 0.18)
DRES3plt <- DRES3 + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Theoretical quantiles of diameter") + 
  scale_x_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3","-2", "-1", "0", "1", "2","3")) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white" ,colour = "black")) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))


DRES4 <- ggplot(fortify(TRDZ, NegativeT))+ stat_qq(size = .7, aes(sample = .stdresid), stroke = .2) + stat_qq_line(aes(sample = .stdresid),size = 0.18) 
DRES4plt <- DRES4 + scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  labs(y = bquote("Effective diffusion standardized residuals"), x = "Theoretical quantiles of zeta-potential") + 
  scale_x_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3), labels = c("-3","-2", "-1", "0", "1", "2","3")) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.4), colour = "black"),
        axis.text.x = element_text( size = rel(.4), colour = "black"),
        axis.ticks = element_line(size = rel(.3)), axis.line = element_blank(),
        axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.4,"cm"),legend.title = element_blank(),
        legend.justification = c(0,0),legend.position = c(.02,.48)) + border("black", size = .65) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white" ,colour = "black"),
        panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) +
  guides(color=guide_legend(override.aes=list(fill=NA))) 


Linearfigurelayout <- egg::ggarrange(RES1plt, RES2plt, RES3plt, RES4plt, SRES1plt, SRES2plt, SRES3plt, SRES4plt, DRES1plt, DRES2plt, DRES3plt, DRES4plt, nrow = 3, ncol = 4)
Linearfigurelayout2 <- egg::ggarrange(SRES1plt,DRES1plt,
                                      SRES2plt,DRES2plt,
                                      SRES3plt,DRES3plt,
                                      SRES4plt,DRES4plt, nrow = 4, ncol = 2)

#ggsave("FIGURESLRanl.pdf", Linearfigurelayout, width = 38, height = 30, units = c("cm"), dpi = 2000, scale = 1:1)
ggsave("Figure_S4_Indepth_SLR_analysis.pdf", Linearfigurelayout2, width = 11.5, height = 21, units = c("cm"), dpi = 2000, scale = 1:1)

SLR_1 = fortify(TRDA, Core_Mucus)

Figure_S4_a_data = data.frame(alpha = SLR_1$alpha, stdresid = SLR_1$.stdresid)

SLR_2 = fortify(TRDDS, SmallMucusT)
Figure_S4_b_data = data.frame(Diameter = SLR_2$Diameter,  stdresid = SLR_2$.stdresid)

SLR_3 = fortify(TRDDA, AlphaMucus1T)
Figure_S4_c_data = data.frame(Diameter = SLR_3$Diameter, stdresid = SLR_3$.stdresid)

SLR_4 = fortify(TRDZ, NegativeT)
Figure_S4_d_data = data.frame(Zeta = SLR_4$Zeta, stdresid = SLR_4$.stdresid)


write.csv(Figure_S4_a_data,'Figure_S4_Indepth_SLR_analysis_a_data.csv', row.names = FALSE)
write.csv(Figure_S4_b_data,'Figure_S4_Indepth_SLR_analysis_b_data.csv', row.names = FALSE)
write.csv(Figure_S4_c_data,'Figure_S4_Indepth_SLR_analysis_c_data.csv', row.names = FALSE)
write.csv(Figure_S4_d_data,'Figure_S4_Indepth_SLR_analysis_d_data.csv', row.names = FALSE)


# Figure_S1_Particle_Type -------------------------------------------------



coohT_ind = which(Mucus_Core$Surface_Chemistry =="COOH")
amineT_ind = which(Mucus_Core$Surface_Chemistry =="AMINE")
pegT_ind = which(Mucus_Core$Surface_Chemistry =="PEG")
chitT_ind = which(Mucus_Core$Surface_Chemistry =="CHITOSAN")
antibodyT_ind = which(Mucus_Core$Surface_Chemistry =="Antibody")
virusT_ind = which(Mucus_Core$Surface_Chemistry == "Virus")
Tmucoad_ind = which(Mucus_Core$Surface_Chemistry == "COOH"     | Mucus_Core$Surface_Chemistry == 'CHITOSAN')
Tmucoin_ind = which(Mucus_Core$Surface_Chemistry == "PEG"      | Mucus_Core$Surface_Chemistry == "AMINE" )
TVA_ind     = which(Mucus_Core$Surface_Chemistry =="Antibody"  | Mucus_Core$Surface_Chemistry == "Virus")

COOHTT       <- Mucus_Core[coohT_ind, ]
AMINETT      <- Mucus_Core[amineT_ind, ]
PEGTT        <- Mucus_Core[pegT_ind, ]
CHITOSANTT   <- Mucus_Core[chitT_ind, ]
ANTIBODIESTT <- Mucus_Core[antibodyT_ind, ]
VirusTT      <- Mucus_Core[virusT_ind, ]
TMucoad      <- Mucus_Core[Tmucoad_ind, ]      
TMucoin      <- Mucus_Core[Tmucoin_ind, ]

Tdiameteron  <- data.frame(Particle = Mucus_Core$Surface_Chemistry, Diameter = Mucus_Core$Diameter)

TChargeon  <- data.frame(Particle = Mucus_Core$Surface_Chemistry, Zeta = Mucus_Core$Zeta)
###
TTempTPon  <- data.frame(Particle = Mucus_Core$Surface_Chemistry, Temperature = Mucus_Core$Temperature)
###
TpHTPon    <- data.frame(Particle = Mucus_Core$Surface_Chemistry, pH = Mucus_Core$pH)
###
TDiffon  <- data.frame(Particle = Mucus_Core$Surface_Chemistry,  Diffusion_constant = Mucus_Core$Diffusion_constant)

Talphaon  <- data.frame(Particle = Mucus_Core$Surface_Chemistry, alpha = Mucus_Core$alpha)
####

###Size
Tboxd <- ggboxplot(Tdiameteron, x = "Particle", y = "Diameter",
                   order = c("Antibody", "Virus", "PEG","AMINE","COOH", "CHITOSAN"),width = .5, lwd = .2)
Tboxpltdia <- Tboxd + scale_y_log10(limits = c(1,1500),breaks = c(1,10,100,1000),  labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Diameter, d (nm)"))+
  border("black", size = .65)+ 
  theme(axis.text = element_text(size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3))) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

###Charge
Tboxz <- ggboxplot(TChargeon, x = "Particle", y = "Zeta",
                   order = c("Antibody", "Virus", "PEG","AMINE", "COOH","CHITOSAN"),width = .5, lwd = .2)
Tboxpltzet <- Tboxz + scale_y_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80"))+ 
  xlab(" ") + ylab(bquote( "Zeta-potential," ~zeta~"(mV)"))+
  border("black", size = .65)+ 
  theme(axis.text = element_text(size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))
### Temperature
TboxTM <- ggboxplot(TTempTPon, x = "Particle", y = "Temperature",
                    order = c("Antibody", "Virus", "PEG","AMINE","COOH", "CHITOSAN"),width = .5, lwd = .2)
TboxpltTem <- TboxTM + scale_y_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310"))+ 
  xlab(" ") + ylab(bquote( "Temperature, T (K)"))+
  border("black", size = .65) + 
  theme(axis.text = element_text(size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))
###pH
TboxpH <- ggboxplot(TpHTPon, x = "Particle", y = "pH",
                    order = c("Antibody", "Virus", "PEG","AMINE","COOH", "CHITOSAN"),width = .5, lwd = .2)
TboxpltpH <- TboxpH + scale_y_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0"))+ 
  xlab(" ") + ylab(bquote( "pH"))+border("black", size = .65) + 
  theme(axis.text = element_text(size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+ 
  theme(axis.ticks = element_line(size = rel(.3)))

###Effective Dif
Tboxdif <- ggboxplot(TDiffon, x = "Particle", y = "Diffusion_constant",
                     order = c("Antibody", "Virus", "PEG","AMINE","COOH", "CHITOSAN"),width = .5, lwd = .2)
Tboxpltdif <- Tboxdif + scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)")) + border("black", size = .65) + 
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3))) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

###alpha
Tboxal <- ggboxplot(Talphaon, x = "Particle", y = "alpha",
                    order = c("Antibody", "Virus", "PEG","AMINE","COOH", "CHITOSAN"),width = .5, lwd = .2)
Tboxpltalp <- Tboxal + scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.00","0.25","0.50","0.75","1.00"))+ 
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+border("black", size = .65) + 
  theme(axis.text = element_text(size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))

TBoxFigurelayout <- egg::ggarrange(Tboxpltdif,Tboxpltalp, 
                                   Tboxpltdia,Tboxpltzet,
                                   TboxpltTem,TboxpltpH, ncol = 2, nrow = 3)
ggsave("Figure_S1_Particle_Type.pdf", TBoxFigurelayout, width = 15, height = 22, units = c("cm"), dpi = 2000, scale = 1:1)

## WRITE TO CSV file
Figure_S1_a_data1 <- data.frame(Antibodies_and_proteins = ANTIBODIESTT$Diffusion_constant)
Figure_S1_a_data2 <- data.frame(Virus                   =      VirusTT$Diffusion_constant)  
Figure_S1_a_data3 <- data.frame(PEG                     =        PEGTT$Diffusion_constant)
Figure_S1_a_data4 <- data.frame(Amine                   =      AMINETT$Diffusion_constant)
Figure_S1_a_data5 <- data.frame(COOH                    =       COOHTT$Diffusion_constant)
Figure_S1_a_data6 <- data.frame(Chitosan                =   CHITOSANTT$Diffusion_constant) 
Figure_S1_a_data <-rbind.fill(Figure_S1_a_data1,Figure_S1_a_data2)
Figure_S1_a_data <-rbind.fill(Figure_S1_a_data,Figure_S1_a_data3)
Figure_S1_a_data <-rbind.fill(Figure_S1_a_data,Figure_S1_a_data4)
Figure_S1_a_data <-rbind.fill(Figure_S1_a_data,Figure_S1_a_data5)
Figure_S1_a_data <-rbind.fill(Figure_S1_a_data,Figure_S1_a_data6)

Figure_S1_b_data1 <- data.frame(Antibodies_and_proteins = ANTIBODIESTT$alpha)
Figure_S1_b_data2 <- data.frame(Virus                   =      VirusTT$alpha)  
Figure_S1_b_data3 <- data.frame(PEG                     =        PEGTT$alpha)
Figure_S1_b_data4 <- data.frame(Amine                   =      AMINETT$alpha)
Figure_S1_b_data5 <- data.frame(COOH                    =       COOHTT$alpha)
Figure_S1_b_data6 <- data.frame(Chitosan                =   CHITOSANTT$alpha) 
Figure_S1_b_data <-rbind.fill(Figure_S1_b_data1,Figure_S1_b_data2)
Figure_S1_b_data <-rbind.fill(Figure_S1_b_data, Figure_S1_b_data3)
Figure_S1_b_data <-rbind.fill(Figure_S1_b_data, Figure_S1_b_data4)
Figure_S1_b_data <-rbind.fill(Figure_S1_b_data, Figure_S1_b_data5)
Figure_S1_b_data <-rbind.fill(Figure_S1_b_data, Figure_S1_b_data6)

Figure_S1_c_data1 <- data.frame(Antibodies_and_proteins = ANTIBODIESTT$Diameter)
Figure_S1_c_data2 <- data.frame(Virus                   =      VirusTT$Diameter)  
Figure_S1_c_data3 <- data.frame(PEG                     =        PEGTT$Diameter)
Figure_S1_c_data4 <- data.frame(Amine                   =      AMINETT$Diameter)
Figure_S1_c_data5 <- data.frame(COOH                    =       COOHTT$Diameter)
Figure_S1_c_data6 <- data.frame(Chitosan                =   CHITOSANTT$Diameter) 
Figure_S1_c_data <-rbind.fill(Figure_S1_c_data1,Figure_S1_c_data2)
Figure_S1_c_data <-rbind.fill(Figure_S1_c_data, Figure_S1_c_data3)
Figure_S1_c_data <-rbind.fill(Figure_S1_c_data, Figure_S1_c_data4)
Figure_S1_c_data <-rbind.fill(Figure_S1_c_data, Figure_S1_c_data5)
Figure_S1_c_data <-rbind.fill(Figure_S1_c_data, Figure_S1_c_data6)

Figure_S1_d_data1 <- data.frame(Antibodies_and_proteins = ANTIBODIESTT$Zeta)
Figure_S1_d_data2 <- data.frame(Virus                   =      VirusTT$Zeta)  
Figure_S1_d_data3 <- data.frame(PEG                     =        PEGTT$Zeta)
Figure_S1_d_data4 <- data.frame(Amine                   =      AMINETT$Zeta)
Figure_S1_d_data5 <- data.frame(COOH                    =       COOHTT$Zeta)
Figure_S1_d_data6 <- data.frame(Chitosan                =   CHITOSANTT$Zeta) 
Figure_S1_d_data <-rbind.fill(Figure_S1_d_data1,Figure_S1_d_data2)
Figure_S1_d_data <-rbind.fill(Figure_S1_d_data, Figure_S1_d_data3)
Figure_S1_d_data <-rbind.fill(Figure_S1_d_data, Figure_S1_d_data4)
Figure_S1_d_data <-rbind.fill(Figure_S1_d_data, Figure_S1_d_data5)
Figure_S1_d_data <-rbind.fill(Figure_S1_d_data, Figure_S1_d_data6)

Figure_S1_e_data1 <- data.frame(Antibodies_and_proteins = ANTIBODIESTT$Temperature)
Figure_S1_e_data2 <- data.frame(Virus                   =      VirusTT$Temperature)  
Figure_S1_e_data3 <- data.frame(PEG                     =        PEGTT$Temperature)
Figure_S1_e_data4 <- data.frame(Amine                   =      AMINETT$Temperature)
Figure_S1_e_data5 <- data.frame(COOH                    =       COOHTT$Temperature)
Figure_S1_e_data6 <- data.frame(Chitosan                =   CHITOSANTT$Temperature) 
Figure_S1_e_data <-rbind.fill(Figure_S1_e_data1,Figure_S1_e_data2)
Figure_S1_e_data <-rbind.fill(Figure_S1_e_data, Figure_S1_e_data3)
Figure_S1_e_data <-rbind.fill(Figure_S1_e_data, Figure_S1_e_data4)
Figure_S1_e_data <-rbind.fill(Figure_S1_e_data, Figure_S1_e_data5)
Figure_S1_e_data <-rbind.fill(Figure_S1_e_data, Figure_S1_e_data6)

Figure_S1_f_data1 <- data.frame(Antibodies_and_proteins = ANTIBODIESTT$pH)
Figure_S1_f_data2 <- data.frame(Virus                   =      VirusTT$pH)  
Figure_S1_f_data3 <- data.frame(PEG                     =        PEGTT$pH)
Figure_S1_f_data4 <- data.frame(Amine                   =      AMINETT$pH)
Figure_S1_f_data5 <- data.frame(COOH                    =       COOHTT$pH)
Figure_S1_f_data6 <- data.frame(Chitosan                =   CHITOSANTT$pH) 
Figure_S1_f_data <-rbind.fill(Figure_S1_f_data1,Figure_S1_f_data2)
Figure_S1_f_data <-rbind.fill(Figure_S1_f_data, Figure_S1_f_data3)
Figure_S1_f_data <-rbind.fill(Figure_S1_f_data, Figure_S1_f_data4)
Figure_S1_f_data <-rbind.fill(Figure_S1_f_data, Figure_S1_f_data5)
Figure_S1_f_data <-rbind.fill(Figure_S1_f_data, Figure_S1_f_data6)


write.csv(Figure_S1_a_data,'Figure_S1_Particle_Type_a_data.csv', row.names = FALSE)
write.csv(Figure_S1_b_data,'Figure_S1_Particle_Type_b_data.csv', row.names = FALSE)
write.csv(Figure_S1_c_data,'Figure_S1_Particle_Type_c_data.csv', row.names = FALSE)
write.csv(Figure_S1_d_data,'Figure_S1_Particle_Type_d_data.csv', row.names = FALSE)
write.csv(Figure_S1_e_data,'Figure_S1_Particle_Type_e_data.csv', row.names = FALSE)
write.csv(Figure_S1_f_data,'Figure_S1_Particle_Type_f_data.csv', row.names = FALSE)
###
mean(alphapoints$D_w)
sd(alphapoints$D_w)

# Mucoadhesive/inert ------------------------------------------------------

Twilmucoaddia <- data.frame(Chemistry = "Mucoadhesive", Diameter           = TMucoad$Diameter)
Twilmucoadzet <- data.frame(Chemistry = "Mucoadhesive", Zeta               = TMucoad$Zeta)
Twilmucoadtemp <- data.frame(Chemistry = "Mucoadhesive", Temperature       = TMucoad$Temperature)
Twilmucoadph <- data.frame(Chemistry = "Mucoadhesive", pH                  = TMucoad$pH)
Twilmucoaddif <- data.frame(Chemistry = "Mucoadhesive", Diffusion_constant = TMucoad$Diffusion_constant)
Twilmucoadalp <- data.frame(Chemistry = "Mucoadhesive", alpha              = TMucoad$alpha)

Twilmucoindia <- data.frame(Chemistry = "Mucoinert",    Diameter           = TMucoin$Diameter)
Twilmucoinzet <- data.frame(Chemistry = "Mucoinert",    Zeta               = TMucoin$Zeta)
Twilmucointemp <- data.frame(Chemistry = "Mucoinert",   Temperature        = TMucoin$Temperature)
Twilmucoinph <- data.frame(Chemistry = "Mucoinert",     pH                 = TMucoin$pH)
Twilmucoindif <- data.frame(Chemistry = "Mucoinert",    Diffusion_constant = TMucoin$Diffusion_constant)
Twilmucoinalp <- data.frame(Chemistry = "Mucoinert",    alpha              = TMucoin$alpha)

Twilmucsdia   <- rbind(Twilmucoaddia, Twilmucoindia)  
Twilmucszet   <- rbind(Twilmucoadzet, Twilmucoinzet)
Twilmucstemp  <- rbind(Twilmucoadtemp, Twilmucointemp)  
TwilmucspH    <- rbind(Twilmucoadph, Twilmucoinph)
Twilmucsdif   <- rbind(Twilmucoaddif, Twilmucoindif) 
Twilmucsalp   <- rbind(Twilmucoadalp, Twilmucoinalp)

wilcox.test(log10(Diameter)~Chemistry,           data =Twilmucsdia)
wilcox.test(Zeta~Chemistry,                      data =Twilmucszet)                       
wilcox.test(log10(Diffusion_constant)~Chemistry, data =Twilmucsdif)
wilcox.test(alpha~Chemistry,                     data =Twilmucsalp)

###
TmucdiaboxT <- ggboxplot(Twilmucsdia, x = "Chemistry", y = "Diameter", 
                         order = c("Mucoadhesive", "Mucoinert"),width = .5, lwd = .2)
Tmucdiaboxplt <- TmucdiaboxT + scale_y_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Diameter, d (nm)"))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
TmuczetboxT <- ggboxplot(Twilmucszet , x = "Chemistry", y = "Zeta", 
                         order = c("Mucoadhesive", "Mucoinert"),width = .5, lwd = .2)
Tmuczetboxplt <- TmuczetboxT + scale_y_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80"))+ 
  xlab(" ") + ylab(bquote( "Zeta-potential,"~zeta~"(mV)"))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
###
TmuctempboxT <- ggboxplot(Twilmucstemp, x = "Chemistry", y = "Temperature", 
                          order = c("Mucoadhesive", "Mucoinert"),width = .5, lwd = .2)
Tmucdtempboxplt <- TmuctempboxT + scale_y_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310"))+ 
  xlab(" ") + ylab(bquote( "Temperature, T (K)"))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+ 
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
TmucpHboxT <- ggboxplot(TwilmucspH , x = "Chemistry", y = "pH", 
                        order = c("Mucoadhesive", "Mucoinert"),width = .5, lwd = .2)
TmucpHboxplt <- TmucpHboxT + scale_y_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0"))+ 
  xlab(" ") + ylab(bquote( "pH"))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))

####
TmucdifboxT <- ggboxplot(Twilmucsdif, x = "Chemistry", y = "Diffusion_constant", 
                         order = c("Mucoadhesive", "Mucoinert"),width = .5, lwd = .2)
Tmucdifboxplt <- TmucdifboxT + scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Effective diffusion,"~D[eff]~"("~mu~m^2~"/s)"))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
TmucalpboxT <- ggboxplot(Twilmucsalp , x = "Chemistry", y = "alpha",
                         order = c("Mucoadhesive", "Mucoinert"),width = .5, lwd = .2)
TmucalpboxplT <- TmucalpboxT + scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.50,.75,1.0), labels = c("0.0", "0.25", "0.50", "0.75","1.0")) +
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  theme(axis.text = element_text(size = rel(.35), color = "Black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))

TfigMUCSTlayoutboxT <- egg::ggarrange(Tmucdifboxplt, TmucalpboxplT,
                                      Tmucdiaboxplt, Tmuczetboxplt,
                                      Tmucdtempboxplt,TmucpHboxplt, 
                                      nrow = 3, ncol = 2)

ggsave("Boxplotadh_iner.pdf", TfigMUCSTlayoutboxT, width = 26, height = 30, units = c("cm"), dpi = 2000)


# Figure_S2_Mucus_Soruce --------------------------------------------------



TMUCTboxDIA <- ggboxplot(Mucus_Core, x = "Mucus_Source", y = "Diameter",
                         order = c("Human cervix","Hydrogel","Pig stomach","Human lung","Pig intestine"),width = .5, lwd = .2)
TMUCTboxpltDIAT <- TMUCTboxDIA + scale_y_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Diameter, d (nm)"))+ border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3))) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

####
TMUCTboxZET <- ggboxplot(Mucus_Core, x = "Mucus_Source", y = "Zeta",order = c("Human cervix","Hydrogel","Pig stomach","Human lung","Pig intestine"),width = .5, lwd = .2)
TMUCTboxpltZET <- TMUCTboxZET  + scale_y_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80"))+ 
  xlab(" ") + ylab(bquote( "Zeta-potential,"~zeta~"(mV)"))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))
###
TMUCTboxTEMP <- ggboxplot(Mucus_Core, x = "Mucus_Source", y = "Temperature",order = c("Human cervix","Hydrogel","Pig stomach","Human lung","Pig intestine"),width = .5, lwd = .2)
TMUCTboxpltTEMP <- TMUCTboxTEMP  + scale_y_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310"))+ 
  xlab(" ") + ylab(bquote( "Temperature, T (K)"))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
TMUCTboxpH <- ggboxplot(Mucus_Core, x = "Mucus_Source", y = "pH",order = c("Human cervix","Hydrogel","Pig stomach","Human lung","Pig intestine"),width = .5, lwd = .2)
TMUCTboxpltpH <- TMUCTboxpH  + scale_y_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0"))+ 
  xlab(" ") + ylab(bquote( "pH"))+ border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
TMUCTboxDIF <- ggboxplot(Mucus_Core, x = "Mucus_Source", y = "Diffusion_constant",order = c("Human cervix","Hydrogel","Pig stomach","Human lung","Pig intestine"),width = .5, lwd = .2)
TMUCTboxpltDIF  <- TMUCTboxDIF+ scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Effective diffusion,"~D[eff]~"("~mu~m^2~"/s)"))+ border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3))) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

####
TMUCTboxALP <- ggboxplot(Mucus_Core, x = "Mucus_Source", y = "alpha",order = c("Human cervix","Hydrogel","Pig stomach","Human lung","Pig intestine"),width = .5, lwd = .2)
TMUCTboxpltALP <- TMUCTboxALP + scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.50,.75,1.0), labels = c("0.0", "0.25", "0.50", "0.75","1.0")) +
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text(size = rel(1), color = "Black"), 
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))

ULTIMATEMUCTYPELAYOUT <- egg::ggarrange(TMUCTboxpltDIF,TMUCTboxpltALP,
                                        TMUCTboxpltDIAT,TMUCTboxpltZET,
                                        TMUCTboxpltTEMP,TMUCTboxpltpH,
                                        nrow = 3,ncol = 2)

ggsave("Figure_S2_Mucus_Source.pdf", ULTIMATEMUCTYPELAYOUT , width = 15, height = 22, units = c("cm"), dpi = 2000)





## Write data into CSV files ##

Figure_S2_a_data1 <- data.frame(Human_cervix    = Mucus_Core[TTISCM_ind,]$Diffusion_constant)
Figure_S2_a_data2 <- data.frame(Hydrogel        = Mucus_Core[TTISMA_ind,]$Diffusion_constant)  
Figure_S2_a_data3 <- data.frame(Pig_stomach     = Mucus_Core[TTISST_ind,]$Diffusion_constant)
Figure_S2_a_data4 <- data.frame(Human_lung      = Mucus_Core[TTISHL_ind,]$Diffusion_constant)
Figure_S2_a_data5 <- data.frame(Pig_intestines  = Mucus_Core[TTISIN_ind,]$Diffusion_constant)

Figure_S2_a_data <-rbind.fill(Figure_S2_a_data1,Figure_S2_a_data2)
Figure_S2_a_data <-rbind.fill(Figure_S2_a_data, Figure_S2_a_data3)
Figure_S2_a_data <-rbind.fill(Figure_S2_a_data, Figure_S2_a_data4)
Figure_S2_a_data <-rbind.fill(Figure_S2_a_data, Figure_S2_a_data5)


Figure_S2_b_data1 <- data.frame(Human_cervix    = Mucus_Core[TTISCM_ind,]$alpha)
Figure_S2_b_data2 <- data.frame(Hydrogel        = Mucus_Core[TTISMA_ind,]$alpha)  
Figure_S2_b_data3 <- data.frame(Pig_stomach     = Mucus_Core[TTISST_ind,]$alpha)
Figure_S2_b_data4 <- data.frame(Human_lung      = Mucus_Core[TTISHL_ind,]$alpha)
Figure_S2_b_data5 <- data.frame(Pig_intestines  = Mucus_Core[TTISIN_ind,]$alpha)

Figure_S2_b_data <-rbind.fill(Figure_S2_b_data1,Figure_S2_b_data2)
Figure_S2_b_data <-rbind.fill(Figure_S2_b_data, Figure_S2_b_data3)
Figure_S2_b_data <-rbind.fill(Figure_S2_b_data, Figure_S2_b_data4)
Figure_S2_b_data <-rbind.fill(Figure_S2_b_data, Figure_S2_b_data5)


Figure_S2_c_data1 <- data.frame(Human_cervix    = Mucus_Core[TTISCM_ind,]$Diameter)
Figure_S2_c_data2 <- data.frame(Hydrogel        = Mucus_Core[TTISMA_ind,]$Diameter)  
Figure_S2_c_data3 <- data.frame(Pig_stomach     = Mucus_Core[TTISST_ind,]$Diameter)
Figure_S2_c_data4 <- data.frame(Human_lung      = Mucus_Core[TTISHL_ind,]$Diameter)
Figure_S2_c_data5 <- data.frame(Pig_intestines  = Mucus_Core[TTISIN_ind,]$Diameter)

Figure_S2_c_data <-rbind.fill(Figure_S2_c_data1,Figure_S2_c_data2)
Figure_S2_c_data <-rbind.fill(Figure_S2_c_data, Figure_S2_c_data3)
Figure_S2_c_data <-rbind.fill(Figure_S2_c_data, Figure_S2_c_data4)
Figure_S2_c_data <-rbind.fill(Figure_S2_c_data, Figure_S2_c_data5)


Figure_S2_d_data1 <- data.frame(Human_cervix    = Mucus_Core[TTISCM_ind,]$Zeta)
Figure_S2_d_data2 <- data.frame(Hydrogel        = Mucus_Core[TTISMA_ind,]$Zeta)  
Figure_S2_d_data3 <- data.frame(Pig_stomach     = Mucus_Core[TTISST_ind,]$Zeta)
Figure_S2_d_data4 <- data.frame(Human_lung      = Mucus_Core[TTISHL_ind,]$Zeta)
Figure_S2_d_data5 <- data.frame(Pig_intestines  = Mucus_Core[TTISIN_ind,]$Zeta)

Figure_S2_d_data <-rbind.fill(Figure_S2_d_data1,Figure_S2_d_data2)
Figure_S2_d_data <-rbind.fill(Figure_S2_d_data, Figure_S2_d_data3)
Figure_S2_d_data <-rbind.fill(Figure_S2_d_data, Figure_S2_d_data4)
Figure_S2_d_data <-rbind.fill(Figure_S2_d_data, Figure_S2_d_data5)


Figure_S2_e_data1 <- data.frame(Human_cervix    = Mucus_Core[TTISCM_ind,]$Temperature)
Figure_S2_e_data2 <- data.frame(Hydrogel        = Mucus_Core[TTISMA_ind,]$Temperature)  
Figure_S2_e_data3 <- data.frame(Pig_stomach     = Mucus_Core[TTISST_ind,]$Temperature)
Figure_S2_e_data4 <- data.frame(Human_lung      = Mucus_Core[TTISHL_ind,]$Temperature)
Figure_S2_e_data5 <- data.frame(Pig_intestines  = Mucus_Core[TTISIN_ind,]$Temperature)

Figure_S2_e_data <-rbind.fill(Figure_S2_e_data1,Figure_S2_e_data2)
Figure_S2_e_data <-rbind.fill(Figure_S2_e_data, Figure_S2_e_data3)
Figure_S2_e_data <-rbind.fill(Figure_S2_e_data, Figure_S2_e_data4)
Figure_S2_e_data <-rbind.fill(Figure_S2_e_data, Figure_S2_e_data5)


Figure_S2_f_data1 <- data.frame(Human_cervix    = Mucus_Core[TTISCM_ind,]$pH)
Figure_S2_f_data2 <- data.frame(Hydrogel        = Mucus_Core[TTISMA_ind,]$pH)  
Figure_S2_f_data3 <- data.frame(Pig_stomach     = Mucus_Core[TTISST_ind,]$pH)
Figure_S2_f_data4 <- data.frame(Human_lung      = Mucus_Core[TTISHL_ind,]$pH)
Figure_S2_f_data5 <- data.frame(Pig_intestines  = Mucus_Core[TTISIN_ind,]$pH)

Figure_S2_f_data <-rbind.fill(Figure_S2_f_data1,Figure_S2_f_data2)
Figure_S2_f_data <-rbind.fill(Figure_S2_f_data, Figure_S2_f_data3)
Figure_S2_f_data <-rbind.fill(Figure_S2_f_data, Figure_S2_f_data4)
Figure_S2_f_data <-rbind.fill(Figure_S2_f_data, Figure_S2_f_data5)



write.csv(Figure_S2_a_data,'Figure_S2_Mucus_source_a_data.csv', row.names = FALSE)
write.csv(Figure_S2_b_data,'Figure_S2_Mucus_source_b_data.csv', row.names = FALSE)
write.csv(Figure_S2_c_data,'Figure_S2_Mucus_source_c_data.csv', row.names = FALSE)
write.csv(Figure_S2_d_data,'Figure_S2_Mucus_source_d_data.csv', row.names = FALSE)
write.csv(Figure_S2_e_data,'Figure_S2_Mucus_source_e_data.csv', row.names = FALSE)
write.csv(Figure_S2_f_data,'Figure_S2_Mucus_source_f_data.csv', row.names = FALSE)


# Purification ------------------------------------------------------------


PURN_ind <-which(Mucus_Core$Purification == "Native")
PURP_ind <-which(Mucus_Core$Purification == "Purified ")
PURP <- Mucus_Core[PURP_ind,]

PURN <- Mucus_Core[PURN_ind,]

PURPDIAon <- data.frame(Purification = PURP$Purification, Diameter = PURP$Diameter)
PURNDIAon <- data.frame(Purification = PURN$Purification, Diameter = PURN$Diameter)
PURDIA    <-  rbind(PURNDIAon,PURPDIAon)

PURPZETon <- data.frame(Purification = PURP$Purification, Zeta = PURP$Zeta)
PURNZETon <- data.frame(Purification = PURN$Purification, Zeta = PURN$Zeta)
PURZET    <-  rbind(PURNZETon,PURPZETon)

PURPTEMon <- data.frame(Purification = PURP$Purification, Temperature = PURP$Temperature)
PURNTEMon <- data.frame(Purification = PURN$Purification, Temperature = PURN$Temperature)
PURTEM <-  rbind(PURNTEMon,PURPTEMon)

PURPpHon <- data.frame(Purification = PURP$Purification, pH = PURP$pH)
PURNpHon <- data.frame(Purification = PURP$Purification, pH = PURP$pH)
PURpH    <-  rbind(PURNpHon,PURPpHon)

PURPDIFon <- data.frame(Purification = PURP$Purification, Diffusion_constant = PURP$Diffusion_constant)
PURNDIFon <- data.frame(Purification = PURN$Purification, Diffusion_constant = PURN$Diffusion_constant)
PURDIF    <-  rbind(PURNDIFon,PURPDIFon)

PURPALPon <- data.frame(Purification = PURP$Purification, alpha = PURP$alpha)
PURNALPon <- data.frame(Purification = PURN$Purification, alpha = PURN$alpha)
PURALP    <-  rbind(PURNALPon,PURPALPon)

PURboxDIA <- ggboxplot(Mucus_Core, x = "Purification", y = "Diameter",width = .5, lwd = .2)
PURboxpltDIA <- PURboxDIA + scale_y_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Diameter, d (nm)"))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
PURboxZET <- ggboxplot(Mucus_Core, x = "Purification", y = "Zeta",width = .5, lwd = .2)
PURboxpltZET <- PURboxZET  + scale_y_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80"))+ 
  xlab(" ") + ylab(bquote( "Zeta-potential,"~zeta~"(mV)"))+theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
###
PURboxTEMP <- ggboxplot(Mucus_Core, x = "Purification", y = "Temperature",width = .5, lwd = .2)
PURboxpltTEMP <- PURboxTEMP  + scale_y_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310"))+ 
  xlab(" ") + ylab(bquote( "Temperature, T (K)"))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
PURboxpH <- ggboxplot(Mucus_Core, x = "Purification", y = "pH",width = .5, lwd = .2)
PURboxpltpH <- PURboxpH  + scale_y_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0"))+ 
  xlab(" ") + ylab(bquote( "pH"))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
PURboxDIF <- ggboxplot(Mucus_Core, x = "Purification", y = "Diffusion_constant",width = .5, lwd = .2)
PURboxpltDIF  <- PURboxDIF+ scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Effective diffusion,"~D[eff]~"("~mu~m^2~"/s)"))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.35), colour = "black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))
####
PURboxALP <- ggboxplot(Mucus_Core, x = "Purification", y = "alpha",width = .5, lwd = .2)
PURboxpltALP <- PURboxALP + scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.50,.75,1.0), labels = c("0.0", "0.25", "0.50", "0.75","1.0")) +
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+ theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45)) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text(size = rel(.35), color = "Black"), 
        axis.title.x = element_text(size = rel(.3)), 
        axis.title.y = element_text(size = rel(.3)))+
  theme(axis.ticks = element_line(size = rel(.3)))

PURfigurelayout <- egg::ggarrange(PURboxpltDIA,PURboxpltZET,PURboxpltTEMP,PURboxpltpH, PURboxpltDIF, PURboxpltALP, nrow = 3, ncol = 2 )
ggsave("BoxplotPurification.pdf", PURfigurelayout , width = 26, height = 30, units = c("cm"), dpi = 2000)


# Figure_S3_Dominant_Mucin ------------------------------------------------




DOMMUC2_ind   <- which(Mucus_Core$Dominate_Mucin == "MUC2")
DOMMUC5AC_ind <- which(Mucus_Core$Dominate_Mucin == "MUC5AC")
DOMMUC5B_ind  <- which(Mucus_Core$Dominate_Mucin == "MUC5B")
DOMMUC2   <- Mucus_Core[DOMMUC2_ind,]
DOMMUC5AC <- Mucus_Core[DOMMUC5AC_ind,]
DOMMUC5B  <- Mucus_Core[DOMMUC5B_ind,]

DOMDIA1 <- data.frame(Dominate_Mucin =   DOMMUC2$Dominate_Mucin, Diameter           =   DOMMUC2$Diameter)
DOMDIA2 <- data.frame(Dominate_Mucin = DOMMUC5AC$Dominate_Mucin, Diameter           = DOMMUC5AC$Diameter)
DOMDIA3 <- data.frame(Dominate_Mucin =  DOMMUC5B$Dominate_Mucin, Diameter           =  DOMMUC5B$Diameter)
DOMDIA  <- rbind(DOMDIA1,DOMDIA2, DOMDIA3)

DOMZET1 <- data.frame(Dominate_Mucin =   DOMMUC2$Dominate_Mucin, Zeta               =   DOMMUC2$Zeta)
DOMZET2 <- data.frame(Dominate_Mucin = DOMMUC5AC$Dominate_Mucin, Zeta               = DOMMUC5AC$Zeta)
DOMZET3 <- data.frame(Dominate_Mucin =  DOMMUC5B$Dominate_Mucin, Zeta               =  DOMMUC5B$Zeta)
DOMZET  <- rbind(DOMZET1,DOMZET2, DOMZET3)

DOMTEM1 <- data.frame(Dominate_Mucin =   DOMMUC2$Dominate_Mucin, Temperature        =   DOMMUC2$Temperature)
DOMTEM2 <- data.frame(Dominate_Mucin = DOMMUC5AC$Dominate_Mucin, Temperature        = DOMMUC5AC$Temperature)
DOMTEM3 <- data.frame(Dominate_Mucin =  DOMMUC5B$Dominate_Mucin, Temperature        =  DOMMUC5B$Temperature)
DOMTEM  <- rbind(DOMTEM1,DOMTEM2, DOMTEM3)

DOMpH1  <- data.frame(Dominate_Mucin =   DOMMUC2$Dominate_Mucin, pH                 =   DOMMUC2$pH)
DOmpH2  <- data.frame(Dominate_Mucin = DOMMUC5AC$Dominate_Mucin, pH                 = DOMMUC5AC$pH)
DOmpH3  <- data.frame(Dominate_Mucin =  DOMMUC5B$Dominate_Mucin, pH                 =  DOMMUC5B$pH)
DOMpH   <- rbind(DOMpH1, DOmpH2, DOmpH3)

DOMDIF1  <- data.frame(Dominate_Mucin =   DOMMUC2$Dominate_Mucin, Diffusion_constant =   DOMMUC2$Diffusion_constant)
DOMDIF2  <- data.frame(Dominate_Mucin = DOMMUC5AC$Dominate_Mucin, Diffusion_constant = DOMMUC5AC$Diffusion_constant)
DOMDIF3  <- data.frame(Dominate_Mucin =  DOMMUC5B$Dominate_Mucin, Diffusion_constant =  DOMMUC5B$Diffusion_constant)
DOMDIF   <- rbind(DOMDIF1,DOMDIF2, DOMDIF3)

DOMALP1 <- data.frame(Dominate_Mucin =   DOMMUC2$Dominate_Mucin, alpha              =   DOMMUC2$alpha)
DOMALP2 <- data.frame(Dominate_Mucin = DOMMUC5AC$Dominate_Mucin, alpha              = DOMMUC5AC$alpha)
DOMALP3 <- data.frame(Dominate_Mucin =  DOMMUC5B$Dominate_Mucin, alpha              =  DOMMUC5B$alpha)
DOMALP  <- rbind(DOMALP1,DOMALP2, DOMALP3) 



DOMboxDIA <- ggboxplot(DOMDIA, x = "Dominate_Mucin", y = "Diameter",
                       order = c("MUC5B","MUC2","MUC5AC"),width = .5, lwd = .2)
DOMboxpltDIA <- DOMboxDIA + scale_y_log10(limits = c(1,1500),breaks = trans_breaks("log10",n= 4, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Diameter, d (nm)"))+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text(size = rel(1), colour = "black"), 
        axis.line = element_blank(),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3))) +
  border("black", size = .65) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

####
DOMboxZET <- ggboxplot(DOMZET , x = "Dominate_Mucin", y = "Zeta",order = c("MUC5B","MUC2","MUC5AC"),width = .5, lwd = .2)
DOMboxpltZET <- DOMboxZET  + scale_y_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80"))+ 
  xlab(" ") + ylab(bquote( "Zeta-potential,"~zeta~"(mV)"))+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.line = element_blank(),
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))+
  border("black", size = .65)
###
DOMboxTEMP <- ggboxplot(DOMTEM, x = "Dominate_Mucin", y = "Temperature",order = c("MUC5B","MUC2","MUC5AC"),width = .5, lwd = .2)
DOMboxpltTEMP <- DOMboxTEMP  + scale_y_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310"))+ 
  xlab(" ") + ylab(bquote( "Temperature, T (K)"))+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.line = element_blank(),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))+
  border("black", size = .65)
####
DOMboxpH <- ggboxplot(DOMpH, x = "Dominate_Mucin", y = "pH",order = c("MUC5B","MUC2","MUC5AC"),width = .5, lwd = .2)
DOMboxpltpH <- DOMboxpH  + scale_y_continuous(limits = c(2.9,8), breaks = c(3,4,5,6,7,8), labels = c("3.0","4.0","5.0","6.0","7.0","8.0"))+ 
  xlab(" ") + ylab(bquote( "pH"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.line = element_blank(),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))+
  border("black", size = .65)
####
DOMboxDIF <- ggboxplot(DOMDIF, x = "Dominate_Mucin", y = "Diffusion_constant",order = c("MUC5B","MUC2","MUC5AC"),width = .5, lwd = .2)
DOMboxpltDIF  <- DOMboxDIF+ scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Effective diffusion,"~D[eff]~"("~mu~m^2~"/s)"))+
  theme(axis.text = element_text( size = rel(1), colour = "black"), 
        axis.line = element_blank(),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))+
  border("black", size = .65) + 
  annotation_logticks(side = "l", size = 0.19, short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +
  # remove clipping
  coord_cartesian(clip="off") 

####
DOMboxALP <- ggboxplot(DOMALP, x = "Dominate_Mucin", y = "alpha",order = c("MUC5B","MUC2","MUC5AC"),width = .5, lwd = .2)
DOMboxpltALP <- DOMboxALP + scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.50,.75,1.0), labels = c("0.0", "0.25", "0.50", "0.75","1.0")) +
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text(size = rel(1), color = "Black"), 
        axis.line = element_blank(),
        axis.text.y = element_text( size = rel(.55), colour = "black"),
        axis.text.x = element_text( size = rel(.55), colour = "black"),
        axis.title.x = element_text(size = rel(.75)), 
        axis.title.y = element_text(size = rel(.75)))+
  theme(axis.ticks = element_line(size = rel(.3)))+
  border("black", size = .65)

DOMfigurelayout <- egg::ggarrange(DOMboxpltDIF, DOMboxpltALP,
                                  DOMboxpltDIA,DOMboxpltZET,
                                  DOMboxpltTEMP,DOMboxpltpH, 
                                  nrow = 3, ncol = 2 )
ggsave("Figure_S3_Dominant_Mucin.pdf", DOMfigurelayout , width = 15, height = 22, units = c("cm"), dpi = 2000)

## Write data into CSV files ##

DDOMMUC2_ind   <- which(DOMDIF$Dominate_Mucin == "MUC2")
DDOMMUC5AC_ind <- which(DOMDIF$Dominate_Mucin == "MUC5AC")
DDOMMUC5B_ind  <- which(DOMDIF$Dominate_Mucin == "MUC5B")

Figure_S3_a_data1 <- data.frame(MUC2          = DOMDIF[DDOMMUC2_ind  ,]$Diffusion_constant)
Figure_S3_a_data2 <- data.frame(MUC5AC        = DOMDIF[DDOMMUC5AC_ind,]$Diffusion_constant)  
Figure_S3_a_data3 <- data.frame(MUC5B         = DOMDIF[DDOMMUC5B_ind ,]$Diffusion_constant)


Figure_S3_a_data <-rbind.fill(Figure_S3_a_data1,Figure_S3_a_data2)
Figure_S3_a_data <-rbind.fill(Figure_S3_a_data, Figure_S3_a_data3)

DDOMMUC2_ind   <- which(DOMALP$Dominate_Mucin == "MUC2")
DDOMMUC5AC_ind <- which(DOMALP$Dominate_Mucin == "MUC5AC")
DDOMMUC5B_ind  <- which(DOMALP$Dominate_Mucin == "MUC5B")


Figure_S3_b_data1 <- data.frame(MUC2          = DOMALP[DDOMMUC2_ind  ,]$alpha)
Figure_S3_b_data2 <- data.frame(MUC5AC        = DOMALP[DDOMMUC5AC_ind,]$alpha)  
Figure_S3_b_data3 <- data.frame(MUC5B         = DOMALP[DDOMMUC5B_ind ,]$alpha)


Figure_S3_b_data <-rbind.fill(Figure_S3_b_data1,Figure_S3_b_data2)
Figure_S3_b_data <-rbind.fill(Figure_S3_b_data, Figure_S3_b_data3)

DDOMMUC2_ind   <- which(DOMDIA$Dominate_Mucin == "MUC2")
DDOMMUC5AC_ind <- which(DOMDIA$Dominate_Mucin == "MUC5AC")
DDOMMUC5B_ind  <- which(DOMDIA$Dominate_Mucin == "MUC5B")

Figure_S3_c_data1 <- data.frame(MUC2           = DOMDIA[DDOMMUC2_ind  ,]$Diameter)
Figure_S3_c_data2 <- data.frame(MUC5AC         = DOMDIA[DDOMMUC5AC_ind,]$Diameter)  
Figure_S3_c_data3 <- data.frame(MUC5B          = DOMDIA[DDOMMUC5B_ind ,]$Diameter)


Figure_S3_c_data <-rbind.fill(Figure_S3_c_data1,Figure_S3_c_data2)
Figure_S3_c_data <-rbind.fill(Figure_S3_c_data, Figure_S3_c_data3)

DDOMMUC2_ind   <- which(DOMZET$Dominate_Mucin == "MUC2")
DDOMMUC5AC_ind <- which(DOMZET$Dominate_Mucin == "MUC5AC")
DDOMMUC5B_ind  <- which(DOMZET$Dominate_Mucin == "MUC5B")

Figure_S3_d_data1 <- data.frame(MUC2           = DOMZET[DDOMMUC2_ind  ,]$Zeta)
Figure_S3_d_data2 <- data.frame(MUC5AC         = DOMZET[DDOMMUC5AC_ind,]$Zeta)  
Figure_S3_d_data3 <- data.frame(MUC5B          = DOMZET[DDOMMUC5B_ind ,]$Zeta)


Figure_S3_d_data <-rbind.fill(Figure_S3_d_data1,Figure_S3_d_data2)
Figure_S3_d_data <-rbind.fill(Figure_S3_d_data, Figure_S3_d_data3)
DDOMMUC2_ind   <- which(DOMTEM$Dominate_Mucin == "MUC2")
DDOMMUC5AC_ind <- which(DOMTEM$Dominate_Mucin == "MUC5AC")
DDOMMUC5B_ind  <- which(DOMTEM$Dominate_Mucin == "MUC5B")


Figure_S3_e_data1 <- data.frame(MUC2          = DOMTEM[DDOMMUC2_ind  ,]$Temperature)
Figure_S3_e_data2 <- data.frame(MUC5AC        = DOMTEM[DDOMMUC5AC_ind,]$Temperature)  
Figure_S3_e_data3 <- data.frame(MUC5B         = DOMTEM[DDOMMUC5B_ind ,]$Temperature)


Figure_S3_e_data <-rbind.fill(Figure_S3_e_data1,Figure_S3_e_data2)
Figure_S3_e_data <-rbind.fill(Figure_S3_e_data, Figure_S3_e_data3)


DDOMMUC2_ind   <- which(DOMpH$Dominate_Mucin == "MUC2")
DDOMMUC5AC_ind <- which(DOMpH$Dominate_Mucin == "MUC5AC")
DDOMMUC5B_ind  <- which(DOMpH$Dominate_Mucin == "MUC5B")

Figure_S3_f_data1 <- data.frame(MUC2           = DOMpH[DDOMMUC2_ind  ,]$pH)
Figure_S3_f_data2 <- data.frame(MUC5AC         = DOMpH[DDOMMUC5AC_ind,]$pH)  
Figure_S3_f_data3 <- data.frame(MUC5B          = DOMpH[DDOMMUC5B_ind ,]$pH)


Figure_S3_f_data <-rbind.fill(Figure_S3_f_data1,Figure_S3_f_data2)
Figure_S3_f_data <-rbind.fill(Figure_S3_f_data, Figure_S3_f_data3)




write.csv(Figure_S3_a_data,'Figure_S3_Dominant_mucin_a_data.csv', row.names = FALSE)
write.csv(Figure_S3_b_data,'Figure_S3_Dominant_mucin_b_data.csv', row.names = FALSE)
write.csv(Figure_S3_c_data,'Figure_S3_Dominant_mucin_c_data.csv', row.names = FALSE)
write.csv(Figure_S3_d_data,'Figure_S3_Dominant_mucin_d_data.csv', row.names = FALSE)
write.csv(Figure_S3_e_data,'Figure_S3_Dominant_mucin_e_data.csv', row.names = FALSE)
write.csv(Figure_S3_f_data,'Figure_S3_Dominant_mucin_f_data.csv', row.names = FALSE)

# 3D Plotting -------------------------------------------------------------


fitDifDIaZet    <- lm(log10(Diffusion_constant)~log10(Diameter)+Zeta, Mucus_Core)
fitAlpDIaZet    <- lm(alpha~log10(Diameter)+Zeta, Mucus_Core)
NfitDifDIaZet   <- lm(log10(Diffusion_constant)~log10(Diameter)+Zeta, NegativeT)
NfitAlpDIaZet   <- lm(alpha~log10(Diameter)+Zeta, NegativeT)
NpHfitDifDIaZet <- lm(log10(Diffusion_constant)~pH+Zeta, NegativeT)
NpHfitAlpDIaZet <- lm(alpha~pH+Zeta, NegativeT)
TFfitDifdiatem  <- lm(log10(Diffusion_constant)~log10(Diameter)+log10(Temperature), Mucus_Core)

summary(fitDifDIaZet)
summary(fitAlpDIaZet)
summary(NfitDifDIaZet)
summary(NfitAlpDIaZet)
summary(NpHfitDifDIaZet)
summary(NpHfitAlpDIaZet)

TfitTDDZ    <- coef(fitDifDIaZet)
TfitTADZ    <- coef(fitAlpDIaZet)
TfitTNDDZ   <- coef(NfitDifDIaZet)
TfitTNADZ   <- coef(NfitAlpDIaZet)
TfitTNpHDDZ   <- coef(NpHfitDifDIaZet)
TfitTNpHADZ   <- coef(NpHfitAlpDIaZet)
TfitTTDDDT <- coef(TFfitDifdiatem)

gridlines = 10.0
xpred  <- seq(min(Mucus_Core$Diameter), max(Mucus_Core$Diameter), by = gridlines)
ypred  <- seq(min(Mucus_Core$Zeta, na.rm = T),  max(Mucus_Core$Zeta, na.rm = T), by = 1)
zpred  <- t(outer(xpred, ypred, function(x,y) 10^(TfitTDDZ[1]+TfitTDDZ[2]*log10(x)+TfitTDDZ[3]*y)))
zapred <- t(outer(xpred, ypred, function(x,y) TfitTADZ[1]+TfitTADZ[2]*log10(x)+TfitTADZ[3]*y))

xpredN  <- seq(min(NegativeT$Diameter),         max(NegativeT$Diameter), by = gridlines)
ypredN  <- seq(min(NegativeT$Zeta, na.rm = T),  max(NegativeT$Zeta, na.rm = T), by = 1)
zpredN  <- t(outer(xpredN, ypredN, function(x,y) 10^(TfitTNDDZ[1]+TfitTNDDZ[2]*log10(x)+TfitTNDDZ[3]*y)))
zapredN <- t(outer(xpredN, ypredN, function(x,y) TfitTNADZ[1]+TfitTNADZ[2]*log10(x)+TfitTNADZ[3]*y))

xpredNpH  <- seq(min(NegativeT$pH, na.rm = T),         max(NegativeT$pH, na.rm = T), by = 0.1)
ypredNpH  <- seq(min(NegativeT$Zeta, na.rm = T),  max(NegativeT$Zeta, na.rm = T), by = 1)
zpredNpH  <- t(outer(xpredNpH, ypredNpH, function(x,y) 10^(TfitTNpHDDZ[1]+TfitTNpHDDZ[2]*x+TfitTNpHDDZ[3]*y)))
zapredNpH <- t(outer(xpredNpH, ypredNpH, function(x,y) TfitTNpHADZ[1]+TfitTNpHADZ[2]*x+TfitTNpHADZ[3]*y))

xpredNTD  <- seq(min(Mucus_Core$Diameter, na.rm = T),         max(Mucus_Core$Diameter, na.rm = T), by = 0.1)
ypredNTD  <- seq(min(Mucus_Core$Temperature, na.rm = T),  max(Mucus_Core$Temperature, na.rm = T), by = 0.1)
zpredNTD  <- t(outer(xpredNTD, ypredNTD, function(x,y) 10^(TfitTTDDDT[1]+TfitTTDDDT[2]*log10(x)+TfitTTDDDT[3]*log10(y))))


###3d plot with plotly
tickproperty <- list(
  size = 14, 
  color = "Black"
)
axisx_prop <- list(
  title = "Diameter",
  tickfont = tickproperty,
  nticks = 2,
  show.grid = T,
  type = "log"
) 
axisxph_prop <- list(
  title = "pH",
  tickfont = tickproperty,
  nticks = 4,
  show.grid = T,
  type = "linear"
) 
axisMC_prop <- list(
  title = "Mucus Concentration",
  tickfont = tickproperty,
  nticks = 4,
  show.grid = T,
  type = "linear"
) 
axisSC_prop <- list(
  title = "Salt Concentration",
  tickfont = tickproperty,
  nticks = 4,
  show.grid = T,
  type = "linear"
) 
axistemp_prop <- list(
  title = "Temperature",
  tickfont = tickproperty,
  nticks = 4,
  show.grid = T,
  type = "log"
) 
axisy_prop <- list(
  title = "Zeta-potential",
  tickfont = tickproperty,
  nticks = 9,
  range = c(-80,80),
  show.grid = T,
  type = "linear"
) 
axisz_prop <- list(
  title = "Effective diffusion",
  tickfont = tickproperty,
  nticks = 5,
  show.grid = T,
  type = "log"
) 
axisza_prop <- list(
  title = "Anomalous exponent",
  tickfont = tickproperty,
  nticks = 4,
  range = c(0,1),
  show.grid = T,
  type = "linear"
) 

scene1 <- list(
  xaxis = axisx_prop,
  yaxis = axisy_prop, 
  zaxis = axisz_prop
)

scene2 <- list(
  xaxis = axisx_prop,
  yaxis = axisy_prop, 
  zaxis = axisza_prop
)

scene3 <- list(
  xaxis = axisxph_prop,
  yaxis = axisy_prop, 
  zaxis = axisz_prop
)

scene4 <- list(
  xaxis = axisxph_prop,
  yaxis = axisy_prop, 
  zaxis = axisza_prop
)
scene5 <- list(
  xaxis = axisx_prop,
  yaxis = axistemp_prop, 
  zaxis = axisz_prop
)
scene6 <- list(
  xaxis = axisy_prop,
  yaxis = axistemp_prop, 
  zaxis = axisz_prop
)
scene7 <- list(
  xaxis = axisSC_prop,
  yaxis = axisMC_prop, 
  zaxis = axisz_prop
)
###Deff
difplot3d <- plot_ly(x = ~xpred, y = ~ypred, z = ~zpred, 
                     type = "surface") %>%
  add_trace(data = Mucus_Core,
            z = ~Diffusion_constant,
            x = ~Diameter,
            y = ~Zeta,
            mode = "markers",
            type = "scatter3d") %>%
  layout(scene = scene1)

difplot3d
##Anomalous exponent
alpplot3d <- plot_ly(x = ~xpred, y = ~ypred, z = ~zapred, 
                     type = "surface") %>%
  add_trace(data = Mucus_Core,
            z = ~alpha,
            x = ~Diameter,
            y = ~Zeta,
            mode = "markers",
            type = "scatter3d") %>%
  layout(scene = scene2)

alpplot3d
##### For ranges that have some statistical significance 

###Deff
Ndifplot3d <- plot_ly(x = ~xpredN, y = ~ypredN, z = ~zpredN, 
                      type = "surface") %>%
  add_trace(data = Mucus_Core,
            z = ~Diffusion_constant,
            x = ~Diameter,
            y = ~Zeta,
            mode = "markers",
            type = "scatter3d") %>%
  layout(scene = scene1)

Ndifplot3d
##Anomalous exponent
Nalpplot3d <- plot_ly(x = ~xpredN, y = ~ypredN, z = ~zapredN, 
                      type = "surface") %>%
  add_trace(data = Mucus_Core,
            z = ~alpha,
            x = ~Diameter,
            y = ~Zeta,
            mode = "markers",
            type = "scatter3d") %>%
  layout(scene = scene2)

Nalpplot3d
##Deff
NpHdifplot3d <- plot_ly(x = ~xpredNpH, y = ~ypredNpH, z = ~zpredNpH, 
                        type = "surface") %>%
  add_trace(data = Mucus_Core,
            z = ~Diffusion_constant,
            x = ~pH,
            y = ~Zeta,
            mode = "markers",
            type = "scatter3d") %>%
  layout(scene = scene3)

NpHdifplot3d
##Anomalous exponent
NpHalpplot3d <- plot_ly(x = ~xpredNpH, y = ~ypredNpH, z = ~zapredNpH, 
                        type = "surface") %>%
  add_trace(data = Mucus_Core,
            z = ~alpha,
            x = ~pH,
            y = ~Zeta,
            mode = "markers",
            type = "scatter3d") %>%
  layout(scene = scene4)

NpHalpplot3d
### Temp and diameter
TDdifplot3d <- plot_ly(x = ~xpredNTD, y= ~ypredNTD, z = ~zpredNTD,
                       type = "surface") %>%
  add_trace(data = Mucus_Core,
            z= ~Diffusion_constant,
            y= ~Temperature,
            x= ~Diameter,
            mode = "markers", 
            type = "scatter3d") %>%
  layout(scene = scene5)
TDdifplot3d

### Temp and zeta
TZdifplot3d <- plot_ly(data = Mucus_Core, x = ~Zeta, y= ~Temperature, z = ~Diffusion_constant, mode = "markers", type = "scatter3d") %>%
  layout(scene = scene6)
TZdifplot3d

##muc and salt conc
SMCdifplot3d <- plot_ly(data = Mucus_Core, x = ~Salt_Concentration, y= ~Muc_Con, z = ~Diffusion_constant, mode = "markers", type = "scatter3d") %>%
  layout(scene = scene7)
SMCdifplot3d


# Random Forest Delete NA -------------------------------------------------

onlyfive <- data.frame(Diffusion_constant = Core_Mucus$Diffusion_constant, alpha = Core_Mucus$alpha, Zeta = Core_Mucus$Zeta, Diameter = Core_Mucus$Diameter,Dominate_Mucin = Core_Mucus$Dominate_Mucin, Mucus_Source = Core_Mucus$Mucus_Source, Surface_Chemistry = Core_Mucus$Surface_Chemistry)


set.seed(888)
Core_Mucus$Mucus_Source <- as.factor(Core_Mucus$Mucus_Source)

str(Core_Mucus)
#RPMCA    <- rfPermute(Diffusion_constant~alpha+Diameter+Zeta+Temperature+pH+Surface_Chemistry+Salt_type+Salt_Concentration+Mucus_Type+          Purification+Dominate_Mucin+MUC2_EL+MUC5AC_EL+MUC5B_EL+MUC6_EL+MUC19_EL,  na.action = na.omit, data = Core_Mucus)
set.seed(888)


#Everything is classified as RF with min. 15 points (alpha,diameter,zeta,pH,Temp.,Surface Chem,Mucus Type,Purification,Dominate Mucin)
set.seed(888)
RPMCA    <- rfPermute(Diffusion_constant~alpha+Zeta+Diameter+Temperature+pH+Surface_Chemistry+Mucus_Source+ Purification+Dominate_Mucin,  na.action = na.omit, data = Core_Mucus)
#Everything without Purification and pH
set.seed(888)
RPMCB    <- rfPermute(Diffusion_constant~alpha+Zeta+Diameter+Temperature+   Surface_Chemistry+Mucus_Source+              Dominate_Mucin,  na.action = na.omit, data = Core_Mucus)
#Important Variables (alpha,zeta,Temperature,surface chemistry,dominate mucin)
set.seed(888)
RPMCC    <- rfPermute(Diffusion_constant~alpha+Zeta+         Mucus_Source+               Surface_Chemistry+Dominate_Mucin                         ,  na.action = na.omit, data = Core_Mucus)
##Combinations of number of elements = 3 between important variables/predictors
set.seed(888)
RPMCD    <- rfPermute(Diffusion_constant~alpha+Zeta+Mucus_Source                                              ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCE    <- rfPermute(Diffusion_constant~alpha+Zeta+Surface_Chemistry                                        ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCF    <- rfPermute(Diffusion_constant~alpha+Zeta+Dominate_Mucin                                           ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCG    <- rfPermute(Diffusion_constant~alpha+     Mucus_Source+Surface_Chemistry                            ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCH    <- rfPermute(Diffusion_constant~alpha+     Mucus_Source+Dominate_Mucin                               ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCI    <- rfPermute(Diffusion_constant~alpha+Surface_Chemistry+Dominate_Mucin                              ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCJ    <- rfPermute(Diffusion_constant~      Zeta+Mucus_Source+Surface_Chemistry                            ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCK    <- rfPermute(Diffusion_constant~      Zeta+Mucus_Source+Dominate_Mucin                               ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCL    <- rfPermute(Diffusion_constant~      Zeta+Surface_Chemistry+Dominate_Mucin                         ,  na.action = na.omit, data = Core_Mucus)
set.seed(888)
RPMCM    <- rfPermute(Diffusion_constant~           Mucus_Source+Surface_Chemistry+Dominate_Mucin             ,  na.action = na.omit, data = Core_Mucus)

### RF hierarchy of Everything
MCTREEe <- rpart(Diffusion_constant~ alpha+Zeta+Diameter+Temperature+pH+Surface_Chemistry+Mucus_Type+Purification+Dominate_Mucin, 
                 method = "anova", data = Core_Mucus, na.action = na.omit,
                 control = rpart.control(minsplit = 10,
                                         maxcompete = 3, cp= 0.001))

fancyRpartPlot(MCTREEe, type = 2)

pruneMCe <- prune(MCTREEe, cp = 0.01)
fancyRpartPlot(pruneMCe)
### RF hierarchy of Everything minus purifaction and pH
MCTREEm <- rpart(Diffusion_constant~ alpha+Zeta+Diameter+Temperature+Surface_Chemistry+Mucus_Type+Dominate_Mucin, 
                 method = "anova", data = Core_Mucus, na.action = na.omit,
                 control = rpart.control(minsplit = 10,
                                         maxcompete = 3, cp= 0.001))

fancyRpartPlot(MCTREEm, type = 2)

pruneMCm <- prune(MCTREEm, cp = 0.01)
fancyRpartPlot(pruneMCm)
### RF hierarchy of Most Important variables
MCTREEi <- rpart(Diffusion_constant~ alpha+Zeta+Mucus_Source+Surface_Chemistry+Dominate_Mucin, 
                 method = "anova", data = Core_Mucus, na.action = na.omit,
                 control = rpart.control(minsplit = 10,
                                         maxcompete = 3, cp= 0.001))

fancyRpartPlot(MCTREEi, type = 2)

pruneMCi <- prune(MCTREEi, cp = 0.01)
fancyRpartPlot(pruneMCi)
###
ImpVar <- c("alpha","Zeta","Temperature","Surface_Chemistry","Dominate Mucin")
ImpVarCom <-combn(ImpVar, 3)

###

RPMCA
rp.importance(RPMCA)
plot(rp.importance(RPMCA), sig.only = FALSE)

RPMCB
rp.importance(RPMCB)
plot(rp.importance(RPMCB), sig.only = FALSE)
plot(rp.importance(RPMCB), n =5, sig.only = FALSE)

RPMCC
rp.importance(RPMCC)
plot(rp.importance(RPMCC), sig.only = FALSE)

RPMCD                       
rp.importance(RPMCD)                       
plot(rp.importance(RPMCD), sig.only = FALSE)

RPMCE
rp.importance(RPMCE)
plot(rp.importance(RPMCE), sig.only = FALSE)

RPMCF
rp.importance(RPMCF)
plot(rp.importance(RPMCF), sig.only = FALSE)

RPMCG
rp.importance(RPMCG)
plot(rp.importance(RPMCG), sig.only = FALSE)

RPMCH
rp.importance(RPMCH)
plot(rp.importance(RPMCH), sig.only = FALSE)

RPMCI
rp.importance(RPMCI)
plot(rp.importance(RPMCI), sig.only = FALSE)

RPMCJ
rp.importance(RPMCJ)
plot(rp.importance(RPMCJ), sig.only = FALSE)

RPMCK
rp.importance(RPMCK)
plot(rp.importance(RPMCK), sig.only = FALSE)

RPMCL
rp.importance(RPMCL)
plot(rp.importance(RPMCL), sig.only = FALSE)

RPMCM
rp.importance(RPMCM)
plot(rp.importance(RPMCM), sig.only = FALSE)

# Random Forest DelNA with log transformations ----------------------------------

set.seed(888)

str(Core_Mucus)
#RPMCA    <- rfPermute(Diffusion_constant~alpha+Diameter+Zeta+Temperature+pH+Surface_Chemistry+Salt_type+Salt_Concentration+Mucus_Type+          Purification+Dominate_Mucin+MUC2_EL+MUC5AC_EL+MUC5B_EL+MUC6_EL+MUC19_EL,  na.action = na.omit, data = Core_Mucus)
set.seed(888)

#Everything is classified as RF with min. 15 points (alpha,diameter,zeta,pH,Temp.,Surface Chem,Mucus Type,Purification,Dominate Mucin)
LPMCA    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Zeta+log10(Core_Mucus$Diameter)+Core_Mucus$Temperature+Core_Mucus$pH+Core_Mucus$Surface_Chemistry+Core_Mucus$Mucus_Source+ Core_Mucus$Purification+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
###Everything minus pH and Purification
set.seed(888)
LPMCB    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Zeta+log10(Core_Mucus$Diameter)+Core_Mucus$Temperature+Core_Mucus$Surface_Chemistry+Core_Mucus$Mucus_Source+                         Core_Mucus$Dominate_Mucin,  na.action = na.omit)
#Important Variables (alpha,zeta,temperature,surface chemistry,dominate mucin)
set.seed(888)
LPMCC    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Zeta+Core_Mucus$Surface_Chemistry+Core_Mucus$Mucus_Source+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
##Combinations of number of elements = 3 between important variables/predictors
set.seed(888)
LPMCD    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Zeta+Core_Mucus$Mucus_Source,  na.action = na.omit)
set.seed(888)
LPMCE    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Zeta+Core_Mucus$Surface_Chemistry,  na.action = na.omit)
set.seed(888)
LPMCF    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Zeta+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
set.seed(888)
LPMCG    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Mucus_Source+Core_Mucus$Surface_Chemistry,  na.action = na.omit)
set.seed(888)
LPMCH    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Mucus_Source+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
set.seed(888)
LPMCI    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$alpha+Core_Mucus$Surface_Chemistry+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
set.seed(888)
LPMCJ    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$Zeta+Core_Mucus$Mucus_Source+Core_Mucus$Surface_Chemistry,  na.action = na.omit)
set.seed(888)
LPMCK    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$Zeta+Core_Mucus$Mucus_Source+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
set.seed(888)
LPMCL    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$Zeta+Core_Mucus$Surface_Chemistry+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
set.seed(888)
LPMCM    <- rfPermute(log10(Core_Mucus$Diffusion_constant)~Core_Mucus$Mucus_Source+Core_Mucus$Surface_Chemistry+Core_Mucus$Dominate_Mucin,  na.action = na.omit)
set.seed(888)
### RF hierarchy of Everything
LCTREEe <- rpart(log10(Diffusion_constant)~ alpha+Zeta+log10(Diameter)+Temperature+pH+Surface_Chemistry+Mucus_Type+Purification+Dominate_Mucin, 
                 method = "anova", data = Core_Mucus, na.action = na.omit,
                 control = rpart.control(minsplit = 10,
                                         maxcompete = 3, cp= 0.001))

fancyRpartPlot(LCTREEe, type = 2)

pruneLCe <- prune(LCTREEe, cp = 0.01)
fancyRpartPlot(pruneLCe)
### RF hierarchy of Everything minus purifaction and pH
LCTREEm <- rpart(log10(Diffusion_constant)~ alpha+Zeta+log10(Diameter)+Temperature+Surface_Chemistry+Mucus_Type+Dominate_Mucin, 
                 method = "anova", data = Core_Mucus, na.action = na.omit,
                 control = rpart.control(minsplit = 10,
                                         maxcompete = 3, cp= 0.001))

fancyRpartPlot(LCTREEm, type = 2)

pruneMCm <- prune(LCTREEm, cp = 0.01)
fancyRpartPlot(pruneMCm)
### RF hierarchy of Most Important variables
LCTREEi <- rpart(log10(Diffusion_constant)~ alpha+Zeta+Mucus_Source+Surface_Chemistry+Dominate_Mucin, 
                 method = "anova", data = Core_Mucus, na.action = na.omit,
                 control = rpart.control(minsplit = 10,
                                         maxcompete = 3, cp= 0.001))

fancyRpartPlot(LCTREEi, type = 2)

pruneMCi <- prune(LCTREEi, cp = 0.01)
fancyRpartPlot(pruneMCi)
###
ImpVar <- c("alpha","Zeta","Mucus_Source","Surface_Chemistry","Dominate Mucin")
ImpVarCom <-combn(ImpVar, 3)

###

LPMCA
rp.importance(LPMCA)
plot(rp.importance(LPMCA), sig.only = FALSE)

LPMCB
rp.importance(LPMCB)
plot(rp.importance(LPMCB), sig.only = FALSE)
plot(rp.importance(LPMCB), n = 5, sig.only = FALSE)

LPMCC
rp.importance(LPMCC)
plot(rp.importance(LPMCC), sig.only = FALSE)

LPMCD                       
rp.importance(LPMCD)                       
plot(rp.importance(LPMCD), sig.only = FALSE)

LPMCE
rp.importance(LPMCE)
plot(rp.importance(LPMCE), sig.only = FALSE)

LPMCF
rp.importance(LPMCF)
plot(rp.importance(LPMCF), sig.only = FALSE)

LPMCG
rp.importance(LPMCG)
plot(rp.importance(LPMCG), sig.only = FALSE)

LPMCH
rp.importance(LPMCH)
plot(rp.importance(LPMCH), sig.only = FALSE)

LPMCI
rp.importance(LPMCI)
plot(rp.importance(LPMCI), sig.only = FALSE)

LPMCJ
rp.importance(LPMCJ)
plot(rp.importance(LPMCJ), sig.only = FALSE)

LPMCK
rp.importance(LPMCK)
plot(rp.importance(LPMCK), sig.only = FALSE)

LPMCL
rp.importance(LPMCL)
plot(rp.importance(LPMCL), sig.only = FALSE)

LPMCM
rp.importance(LPMCM)
plot(rp.importance(LPMCM), sig.only = FALSE)



# Figure_1_Random_forest --------------------------------------------------

# % IncMSE 
RFimpvec <- c(22,19,13,15,10)
# Groups in order with respective average %IncMSE
RFimpvec2 <- c("Anomalous exponent, ","Particle type","Mucus source","Zeta potential, ","Dominant mucin" )
# Standard error with each respective group
RFimpvecstd <- c(3,7,7,5,8)
# p-value of each variable
RFimpvec3 <- c(0.01,0.02,0.05,0.1,0.1)
#Combination of data frame
RFimpplot <- data.frame(INCMSE = RFimpvec, Importance = RFimpvec2, Std = RFimpvecstd, pval = RFimpvec3)


RFp1 <- ggplot(RFimpplot, aes(x = reorder(Importance, -INCMSE), y = INCMSE))
RFplt1 <- RFp1  + geom_bar(stat = "identity", width = 0.5, fill = "gray") + 
  geom_errorbar(aes(ymin = INCMSE-Std, ymax = INCMSE+Std), width = 0.1, size = .2)+
  xlab( " ") + ylab("% MSE Increase")+
  scale_y_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,25,30), labels = c("0","5","10","15","20","25","30"))+
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, vjust = .71))+
  theme(axis.ticks = element_line(size = rel(.3), color = "Black"))+ 
  theme(axis.title.y = element_text(size = rel(.45)))+
  theme(axis.line = element_blank())+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+  
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.border = element_rect(color = "black"))+
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))
#egg::ggarrange(RFplt1, NULL,nrow = 1, ncol = 2)
RFfigurelayout <- plot_grid(RFplt1, NULL, ncol = 2)
ggsave("Figure_1_RFMSE.pdf", RFfigurelayout, width = 13, height = 7, units = c("cm"), dpi = 2000, scale = 1:1)

pdf("Figure_1_rpart.pdf", width = 13/2.54, height = 7/2.54,useDingbats=FALSE)
par(mfrow = c(1,2))
par(mfg = c(1,1))
empty.plot()
rpart.plot(LCTREEi,
           type = 4, 
           extra = 100,
           clip.right.labs = FALSE,
           under = TRUE,
           prefix = "log Deff um^2/s\n h\n", 
           tweak = .95,
           under.cex = 1.05)
dev.off()
par(mfrow = c(1,1))


grid.arrange(rpart.plot(LCTREEi),RFplt1)


# RF missForest -----------------------------------------------------------
MFore_Mucus <- Core_Mucus
MFore_Mucus[,4] <- NULL 
MFore_Mucus[,5] <- NULL
MFore_Mucus[,23] <- NULL
MRFC <- missForest(MFore_Mucus, ntree = 500, verbose = TRUE)
CFRM <- MRFC$ximp

##Everything
set.seed(888)
CMPRA <- rfPermute(Diffusion_constant~alpha+Zeta+Diameter+Temperature+pH+Surface_Chemistry+Salt_type+Salt_Concentration+Mucus_Type+Muc_Con+Purification+Dominate_Mucin, na.action = na.omit, data = CFRM)
set.seed(888)
LMPRA <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+CFRM$Zeta+log10(CFRM$Diameter)+log10(CFRM$Temperature)+CFRM$pH+CFRM$Surface_Chemistry+CFRM$Salt_type+CFRM$Salt_Concentration+CFRM$Mucus_Type+CFRM$Muc_Con+CFRM$Purification+CFRM$Dominate_Mucin, na.action = na.omit)
### Everything minus non sign MSE pval for both log and non-log transformation (minus pH, Muc_Con, Purification)
set.seed(888)
CMPRB <- rfPermute(Diffusion_constant~alpha+Zeta+Diameter+Temperature+   Surface_Chemistry+Salt_type+Salt_Concentration+Mucus_Type+                     Dominate_Mucin, na.action = na.omit, data = CFRM)
set.seed(888)
LMPRB <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+CFRM$Zeta+log10(CFRM$Diameter)+log10(CFRM$Temperature)+        CFRM$Surface_Chemistry+CFRM$Salt_type+CFRM$Salt_Concentration+CFRM$Mucus_Type+                               CFRM$Dominate_Mucin, na.action = na.omit)
#### 5 Most Important variable based on their own permute
### For no log transforamtion (alpha,diameter,surface chem, mucus type and dominate mucin)
set.seed(888)
CMPRC <- rfPermute(Diffusion_constant~alpha+Diameter+   Surface_Chemistry+Mucus_Type+                     Dominate_Mucin, na.action = na.omit, data = CFRM)
set.seed(888)
Ccom <- c("alpha","diameter","surface Chemistry","Mucus Type", "Dominate Mucin")
set.seed(888)
CcomC <- combn(Ccom, 3)
set.seed(888)
CMPRD <- rfPermute(Diffusion_constant~alpha+Diameter+Surface_Chemistry,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRE <- rfPermute(Diffusion_constant~alpha+Diameter+Mucus_Type,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRF <- rfPermute(Diffusion_constant~alpha+Diameter+Dominate_Mucin,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRG <- rfPermute(Diffusion_constant~alpha+Surface_Chemistry+Mucus_Type,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRH <- rfPermute(Diffusion_constant~alpha+Surface_Chemistry+Dominate_Mucin,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRI <- rfPermute(Diffusion_constant~alpha+Mucus_Type+Dominate_Mucin,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRJ <- rfPermute(Diffusion_constant~Diameter+Surface_Chemistry+Mucus_Type,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRK <- rfPermute(Diffusion_constant~Diameter+Surface_Chemistry+Dominate_Mucin,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRL <- rfPermute(Diffusion_constant~Diameter+Mucus_Type+Dominate_Mucin,na.action = na.omit, data = CFRM)
set.seed(888)
CMPRM <- rfPermute(Diffusion_constant~Surface_Chemistry+Mucus_Type+Dominate_Mucin,na.action = na.omit, data = CFRM)

### For log traFnsformation (alpha, zeta, diameter,surface chem, mucus type)
LMPRC <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+CFRM$Zeta+log10(CFRM$Diameter)+        CFRM$Surface_Chemistry+CFRM$Mucus_Type                               , na.action = na.omit)
Lcom <- c("alpha","zeta","diameter","Surface Chemistry", "Mucus Type")
LcomC <- combn(Lcom, 3)
set.seed(888)
LMPRD    <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+CFRM$Zeta+log10(CFRM$Diameter)                                                                                ,  na.action = na.omit)
set.seed(888)
LMPRE    <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+CFRM$Zeta+                           CFRM$Surface_Chemistry                                                   ,  na.action = na.omit)
set.seed(888)
LMPRF    <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+CFRM$Zeta+                                                        CFRM$Mucus_Type                             ,  na.action = na.omit)
set.seed(888)
LMPRG    <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+                log10(CFRM$Diameter)+CFRM$Surface_Chemistry                                                   ,  na.action = na.omit)
set.seed(888)
LMPRH    <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+                log10(CFRM$Diameter)+                             CFRM$Mucus_Type                             ,  na.action = na.omit)
set.seed(888)
LMPRI    <- rfPermute(log10(CFRM$Diffusion_constant)~CFRM$alpha+                                           CFRM$Surface_Chemistry+CFRM$Mucus_Type                             ,  na.action = na.omit)
set.seed(888)
LMPRJ    <- rfPermute(log10(CFRM$Diffusion_constant)~                 CFRM$Zeta+log10(CFRM$Diameter)+CFRM$Surface_Chemistry                                                   ,  na.action = na.omit)
set.seed(888)
LMPRK    <- rfPermute(log10(CFRM$Diffusion_constant)~                 CFRM$Zeta+log10(CFRM$Diameter)+                             CFRM$Mucus_Type                             ,  na.action = na.omit)
set.seed(888)
LMPRL    <- rfPermute(log10(CFRM$Diffusion_constant)~                 CFRM$Zeta+                           CFRM$Surface_Chemistry+CFRM$Mucus_Type                             ,  na.action = na.omit)
set.seed(888)
LMPRM    <- rfPermute(log10(CFRM$Diffusion_constant)~                                 log10(CFRM$Diameter)+CFRM$Surface_Chemistry+CFRM$Mucus_Type                             ,  na.action = na.omit)


CMPRA
CMPRB
CMPRC
CMPRD
CMPRE
CMPRF
CMPRG
CMPRH
CMPRI
CMPRJ
CMPRK
CMPRL
CMPRM

rp.importance(CMPRA)
rp.importance(CMPRB)
rp.importance(CMPRC)
rp.importance(CMPRD)
rp.importance(CMPRE)
rp.importance(CMPRF)
rp.importance(CMPRG)
rp.importance(CMPRH)
rp.importance(CMPRI)
rp.importance(CMPRJ)
rp.importance(CMPRK)
rp.importance(CMPRL)
rp.importance(CMPRM)




plot(rp.importance(CMPRA), sig.only = FALSE)
plot(rp.importance(CMPRB), sig.only = FALSE)
plot(rp.importance(CMPRB), n = 5, sig.only = FALSE)
plot(rp.importance(CMPRC), sig.only = FALSE)
plot(rp.importance(CMPRD), sig.only = FALSE)
plot(rp.importance(CMPRE), sig.only = FALSE)
plot(rp.importance(CMPRF), sig.only = FALSE)
plot(rp.importance(CMPRG), sig.only = FALSE)
plot(rp.importance(CMPRH), sig.only = FALSE)
plot(rp.importance(CMPRI), sig.only = FALSE)
plot(rp.importance(CMPRJ), sig.only = FALSE)
plot(rp.importance(CMPRK), sig.only = FALSE)
plot(rp.importance(CMPRL), sig.only = FALSE)
plot(rp.importance(CMPRM), sig.only = FALSE)


### log transformations

LMPRA
LMPRB
LMPRC
LMPRD
LMPRE
LMPRF
LMPRG
LMPRH
LMPRI
LMPRJ
LMPRK
LMPRL
LMPRM

rp.importance(LMPRA)
rp.importance(LMPRB)
rp.importance(LMPRC)
rp.importance(LMPRD)
rp.importance(LMPRE)
rp.importance(LMPRF)
rp.importance(LMPRG)
rp.importance(LMPRH)
rp.importance(LMPRI)
rp.importance(LMPRJ)
rp.importance(LMPRK)
rp.importance(LMPRL)
rp.importance(LMPRM)


plot(rp.importance(LMPRA), sig.only = FALSE)
plot(rp.importance(LMPRB), sig.only = FALSE)
plot(rp.importance(LMPRB), n = 5, sig.only = FALSE)
plot(rp.importance(LMPRC), sig.only = FALSE)
plot(rp.importance(LMPRD), sig.only = FALSE)
plot(rp.importance(LMPRE), sig.only = FALSE)
plot(rp.importance(LMPRF), sig.only = FALSE)
plot(rp.importance(LMPRG), sig.only = FALSE)
plot(rp.importance(LMPRH), sig.only = FALSE)
plot(rp.importance(LMPRI), sig.only = FALSE)
plot(rp.importance(LMPRJ), sig.only = FALSE)
plot(rp.importance(LMPRK), sig.only = FALSE)
plot(rp.importance(LMPRL), sig.only = FALSE)
plot(rp.importance(LMPRM), sig.only = FALSE)


# Plots based on missForest Data ------------------------------------------
dif_ind <- which(CFRM$Diffusion_constant < 6)
dif_ind2 <- which(CFRM$Diffusion_constant > 6)

missCFRM <- CFRM[dif_ind,] #non constant alpha values
missCFRM2 <- CFRM[dif_ind2,] #constant alpha values

RFm1 <- ggplot(CFRM, aes(x= alpha, y=Diffusion_constant))+geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
RFm1plt <- RFm1 + scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) + scale_y_log10(limits = c(0.1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"), x = "Anomalous exponent,"~alpha) + 
  scale_x_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"), axis.ticks = element_line(size = rel(.3)),
        axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4)))+
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.3,"line"),
        legend.spacing = unit(0.1,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.08)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_smooth(data = missCFRM, method = "lm", se = T, show.legend = F, level = 0.95, color = "black", size = .18, na.rm = TRUE) +
  geom_smooth(data = missCFRM2, method = "lm", se = F, show.legend = F, level = 0.95, color = "black", size = .18, na.rm = TRUE) +
  border("black", size = .65) +
  coord_cartesian(xlim = c(0,1.1), ylim = c(1e-5,1000))

RFt1 <- lm(log10(Diffusion_constant)~alpha, data = missCFRM)
summary(RFt1)

cor.test(Core_Mucus$alpha,log10(CFRM$Diffusion_constant), 
         method = "spearman",use = "na.or.complete", exact = FALSE)
cor.test(Core_Mucus$alpha,log10(CFRM$Diffusion_constant), 
         method = "pearson",use = "na.or.complete", exact = FALSE)


missalb <- lm(alpha~log10(Diffusion_constant), data = missCFRM)
summary(missalb)

RFb <- ggplot(data = CFRM, aes(x = Diffusion_constant, y = alpha)) + geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)
RFbplt <- RFb + scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies 
and proteins","Chitosan","COOH", "PEG","Viruses")) + scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) +
  scale_x_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  theme(axis.text = element_text( size = rel(.4), colour = "black"), axis.ticks = element_line(size = rel(.3)),
        axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4)))+
  labs(x = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"), y = "Anomalous exponent,"~alpha) +
  theme(legend.key = element_blank(), legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"), legend.key.width = unit(.2,"line"),
        legend.spacing = unit(0.1,"cm"),legend.title = element_blank(), legend.justification = c(0,0),legend.position = c(.02,.48)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white" ,colour = "black")) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_smooth(data = missCFRM, method = "lm", se = T, show.legend = F, level = 0.95, color = "black", size = .18, na.rm = TRUE) +
  geom_smooth(data = missCFRM2, method = "lm", se = F, show.legend = F, level = 0.95, color = "black", size = .18, na.rm = TRUE) + border("black", size = .65) 
#geom_vline(xintercept = 3.6, linetype = "dotted", size = 0.2)

RFafiglayout <- egg::ggarrange(RFm1plt, RFbplt, nrow = 1)

ggsave("missForestDeffvalpha.pdf", RFafiglayout, width = 10, height = 5, units = c("cm"), dpi = 2000)



RFTc1 <- ggplot(CFRM, aes(x = Zeta, y = Diffusion_constant))

RFTc1plt <- RFTc1 + geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2) +
  scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies              
 and proteins","Chitosan","COOH", "PEG","Viruses"))+
  scale_x_continuous(limits = c(-80,80), breaks = c(-80,-60,-40,-20,0,20,40,60,80), labels = c("-80","-60","-40","-20","0","20","40","60","80")) +
  scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))  +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+ theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4))) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Zeta-potential,"~zeta~"(mV)")+theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank())+ theme(legend.text = element_text(size = 3.65),legend.background = element_rect(fill = "transparent"))+ theme(legend.key.width = unit(.3,"line"))+
  theme(legend.spacing = unit(0.4,"cm")) + theme(legend.title = element_blank()) + theme(legend.justification = c(0,0),legend.position = c(.02,.48))+ theme(panel.grid.major = element_blank())+
  geom_smooth(data = CFRM[which(CFRM$Zeta < 0),], method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid") +
  guides(color=guide_legend(override.aes=list(fill=NA)))+theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = "dotted", size =.2)+ border("black", size = .65)

RFTc2 <- ggplot(CFRM, aes(x = Zeta, y = alpha))+geom_point(size = .7, aes(shape = Surface_Chemistry),stroke = .2)+ scale_shape_manual(values = c(6,4,2,1,0,5))
RFTc2plt <- RFTc2  + 
  theme(axis.title.x = element_text(size = rel(.4)), axis.title.y = element_text(size = rel(.4)) ) + theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha), x = bquote("Zeta-potential,"~zeta~"(mV)")) + theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80,-60,-40,-20, 0, 20, 40, 60, 80) ) + geom_vline(xintercept = 0, linetype = "dotted", size =.2) + 
  scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) +theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(axis.ticks = element_line(size = rel(.3))) + border("black", size = .65) +
  geom_smooth(data = NegativeT, method = "lm", se = T, show.legend = F, na.rm = TRUE, color = "black", size = .18, linetype = "solid")

Tfiglayout2 <- egg::ggarrange(RFTc1plt, RFTc2plt, nrow = 1)

ggsave("missForestCT.pdf", Tfiglayout2, width = 10, height = 5, units = c("cm"), dpi = 2000, scale = 1:1)


# multivariable analysis --------------------------------------------------

MVA1 <- lm(log10(Diffusion_constant)~log10(Diameter)+Zeta, Mucus_Core)
MVA2 <- lm(log10(Diffusion_constant)~log10(Diameter)+log10(Temperature), Mucus_Core)
MVA3 <- lm(log10(Diffusion_constant)~Zeta+log10(Temperature), Mucus_Core)
MVA4 <- lm(log10(Diffusion_constant)~log10(Diameter)+ Muc_Con, Mucus_Core)
MVA5 <- lm(log10(Diffusion_constant)~log10(Diameter)+Salt_Concentration, Mucus_Core)
MVA6 <- lm(log10(Diffusion_constant)~Zeta + Muc_Con, Mucus_Core)
MVA7 <- lm(log10(Diffusion_constant)~Zeta + Salt_Concentration, Mucus_Core)
MVA8 <- lm(log10(Diffusion_constant)~pH+ Muc_Con, Mucus_Core)
MVA9 <- lm(log10(Diffusion_constant)~pH+Salt_Concentration, Mucus_Core)

MVAN1 <- lm(log10(Diffusion_constant)~log10(Diameter)+Zeta, NegativeT)
MVAN2 <- lm(log10(Diffusion_constant)~pH+Zeta, NegativeT)



MVAa1 <- lm(alpha~log10(Diameter)+Zeta, Mucus_Core)
MVAa2 <- lm(alpha~log10(Diameter)+Zeta, NegativeT)
MVAa3 <- lm(alpha~pH+Zeta, NegativeT)

summary(MVA1) 
summary(MVA2) 
summary(MVA3) 
summary(MVA4) 
summary(MVA5) 
summary(MVA6) 
summary(MVA7) 
summary(MVA8) 
summary(MVA9) 
summary(MVAN1) 
summary(MVAN2) 
summary(MVAa1)
summary(MVAa2)
summary(MVAa3)


# Temporary Temperature graph placement -----------------------------------


TT1 <- ggplot(Mucus_Core, aes(x = Temperature, y = Diffusion_constant))+
  geom_point(size = .7, aes(shape = Surface_Chemistry), stroke = .2)

TT1plt <- TT1  +
  scale_shape_manual(values = c(6,4,2,1,0,5), labels = c("AMINE","Antibodies              
                                                         and proteins","Chitosan","COOH", "PEG","Viruses"))+
  scale_x_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310")) +
  scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))  +
  theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  theme(axis.ticks = element_line(size = rel(.3)))+
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45))) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Effective diffusion," ~D[eff]~"("~mu~m^2~"/s)"),x = "Temperature, T (K)")+
  theme(legend.key.height = unit(.5, "line"))+
  theme(legend.key = element_blank())+
  theme(legend.text = element_text(size =  rel(.45)),legend.background = element_rect(fill = "transparent"))+ 
  theme(legend.key.width = unit(.3,"line"))+
  theme(legend.spacing = unit(0.4,"cm")) + 
  theme(legend.title = element_blank()) + theme(legend.justification = c(0,0),legend.position = c(.02,.48))+ 
  theme(panel.grid.major = element_blank())+ 
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.grid.minor = element_blank()) 

TT2 <- ggplot(Mucus_Core, aes(x = Temperature, y = alpha))+geom_point(size = .7, aes(shape = Surface_Chemistry, color = Data_type),stroke = .2)+ scale_shape_manual(values = c(6,4,2,1,0,5))
TT2plt <- TT2 + scale_color_manual(values = c("Black", "lightsteelblue3"), labels = c("Experiment", "Prediction")) +
  theme(axis.title.x = element_text(size = rel(.45)), axis.title.y = element_text(size = rel(.45)) ) + 
  theme(panel.background = element_rect(fill = "white" ,colour = "black"))+
  labs(y = bquote("Anomalous exponent,"~alpha), x = bquote("Temperature, T (K)")) + 
  theme(legend.position = "none") +theme(axis.text = element_text(size = rel(.4), colour = "black"))+
  scale_x_continuous(limits = c(294,311), breaks = c(295,300,305,310),  labels =c("295","300","305","310")) +  
  scale_y_continuous(limits = c(0,1.1), breaks = c(0.0,0.25,0.50,0.75,1.0), labels = c("0.0","0.25","0.50","0.75","1.0")) +
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(colour = "black", linetype = "solid", size = 0.45))+
  theme(axis.ticks = element_line(size = rel(.3))) 

TTEMPboxDIF <- ggboxplot(TempTboxD, x = "Temp_group", y = "Diffusion_constant",width = .5, lwd = .2)
TTEMPboxpltDIF  <- TTEMPboxDIF+ scale_y_log10(limits = c(1E-5, 1000), breaks = trans_breaks("log10",n= 8, function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+ 
  xlab(" ") + ylab(bquote( "Effective diffusion,"~D[eff]~"("~mu~m^2~"/s)"))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text( size = rel(.4), colour = "black"), 
        axis.title.x = element_text(size = rel(.45)), 
        axis.title.y = element_text(size = rel(.45))) + 
  theme(axis.ticks = element_line(size = rel(.3)))
####
TTEMPboxALP <- ggboxplot(TempTboxA, x = "Temp_group", y = "alpha", width = .5, lwd = .2) 
TTEMPboxpltALP <- TTEMPboxALP + scale_y_continuous(limits = c(0,1.1), breaks = c(0,.25,.50,.75,1.0), labels = c("0.0", "0.25", "0.50", "0.75","1.0")) +
  xlab(" ") + ylab(bquote("Anomalous exponent,"~alpha))+border("black", size = .65) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text = element_text(size = rel(.4), color = "Black"), 
        axis.title.x = element_text(size = rel(.45)), 
        axis.title.y = element_text(size = rel(.45))) + 
  theme(axis.ticks = element_line(size = rel(.3)))


TfiglayoutT <- egg::ggarrange(TT1plt, TT2plt,TTEMPboxpltDIF,TTEMPboxpltALP, nrow = 2, ncol = 2)
ggsave("FIGURETT.pdf", TfiglayoutT, width = 14, height = 12, units = c("cm"), dpi = 2000, scale = 1:1)



# END OF CODE -------------------------------------------------------------
sink()


closeAllConnections()
