rm(list=ls())
# --------------------------------------------------------------------------------------------------------
# (1) Install and load required packages -----------------------------------------------------------------
if (!'lme4' %in% installed.packages()) install.packages('lme4') ; require(lme4)
if (!'car' %in% installed.packages()) install.packages('car') ; require(car)
if (!'MuMIn' %in% installed.packages()) install.packages('MuMIn') ; require(MuMIn)
if (!'tidyverse' %in% installed.packages()) install.packages('tidyverse') ; require(tidyverse)
if (!'plyr' %in% installed.packages()) install.packages('plyr') ; require(plyr)
if (!'broom' %in% installed.packages()) install.packages('broom') ; require(broom)
if (!'coda' %in% installed.packages()) install.packages('coda') ; require(coda)
if (!'grid' %in% installed.packages()) install.packages('grid') ; require(grid)
if (!'gridExtra' %in% installed.packages()) install.packages('gridExtra') ; require(gridExtra)
if (!'brms' %in% installed.packages()) install.packages('brms') ; require(brms)
if (!'broom.mixed' %in% installed.packages()) install.packages('broom.mixed') ; require(broom.mixed)
if (!'merTools' %in% installed.packages()) install.packages('merTools') ; require(merTools)
if (!'tidybayes' %in% installed.packages()) install.packages('tidybayes') ; require(tidybayes)
if (!'parallel' %in% installed.packages()) install.packages('parallel') ; require(parallel)
if (!'rptR' %in% installed.packages()) install.packages('rptR') ; require(rptR)
if (!'RColorBrewer' %in% installed.packages()) install.packages('RColorBrewer') ; require(RColorBrewer)
display.brewer.all(n=4, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=TRUE)
display.brewer.pal(10, "RdBu")

col.born.capt <- brewer.pal(n = 10, name = "RdBu")[3]
col.born.wild <- brewer.pal(n = 10, name = "RdBu")[2]
col.suma.capt <- brewer.pal(n = 10, name = "RdBu")[8]
col.suma.wild <- brewer.pal(n = 10, name = "RdBu")[9]
 
# --------------------------------------------------------------------------------------------------------
# (2) Read in all data & settings ------------------------------------------------------------------------
path <- "D:\\R\\MS Orang IndDiff\\"
dat <- read.table(paste0(path,"data_plastic_orang.csv"), header=TRUE, sep=",", stringsAsFactors=TRUE)

contr <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))

# extract legend function
g_legend <- function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}

# --------------------------------------------------------------------------------------------------------
# (3) Prepare data ---------------------------------------------------------------------------------------
# Subset data for 1st variable: Persistence (mothers' signals versus responses) ---
test.data <- as.data.frame(na.omit(dat[ ,c("Persist", "Context","parity", "inf_age","kinship","agediff", "species", "specset", "group", "setting","ID_obs", "animal_id", "ID_rec")]))
test.data$context_jt <- as.numeric(test.data$Context==levels(test.data$Context)[3])
test.data.sig <- subset(test.data, kinship=='motoff' & agediff =='younger')
test.data.sig$context_jt <- as.factor(test.data.sig$context_jt)

# Omit all levels of animal_id that contributed fewer than 30 cases
NID <- data.frame(table(test.data.sig$animal_id)) ; colnames(NID) <- c("animal_id","N")
NID$animal_id <- as.character(NID$animal_id)
NID30 <- subset(NID, N>=2)
test.data.sig <- subset(test.data.sig, animal_id %in% NID30$animal_id)

# Subset data for 2nd variable: Responsiveness to infant requests --- 
test.data2 <- as.data.frame(na.omit(dat[ ,c("ASO", "Context","parity", "inf_age","kinship","agediff", "species","specset", "group", "setting","ID_obs", "animal_id", "ID_rec")]))
test.data2$context_fs <- as.numeric(test.data2$Context==levels(test.data2$Context)[1])
test.data.out <- subset(test.data2, kinship=='motoff' & agediff =='older')

# Omit all levels of animal_id that contributed fewer than 30 cases
NID <- data.frame(table(test.data.out$ID_rec)) ; colnames(NID) <- c("ID_rec","N")
NID$ID_rec <- as.character(NID$ID_rec)
NID30 <- subset(NID, N>=2)
test.data.out <- subset(test.data.out, ID_rec %in% NID30$ID_rec)

# --------------------------------------------------------------------------------------------------------
# (4) Check collinearity ---------------------------------------------------------------------------------
vif(lm(rnorm(nrow(test.data.sig)) ~ setting + species + inf_age + parity + context_jt, data = test.data.sig))
vif(lm(rnorm(nrow(test.data.out)) ~ setting + species + inf_age + parity + context_fs, data = test.data.out))
### collinearity check ok: max vif = 1.7

# --------------------------------------------------------------------------------------------------------
# (5) Analyses: 1st variable: Persistence ----------------------------------------------------------------
# ----------------------------------------------------
# (5.1) Frequentist approach -------------------------
# Random intercept model
ri <- glmer(formula = Persist ~ setting + species + inf_age + parity + context_jt +
	(1|animal_id) + (1|group), family = binomial, data = test.data.sig, control = contr)

# Random slope model
ris <- glmer(formula = Persist ~ setting + species + inf_age + parity + context_jt +
	(1|animal_id) + (1|group) + (0+context_jt|animal_id), family = binomial, data = test.data.sig, control = contr)

length(residuals(ri)) # 650
length(residuals(ris)) # 650

# Check models with lmer 
rin <- lmer(formula = Persist ~ setting + species + inf_age + parity + context_jt +
	(1|animal_id) + (1|group), data = test.data.sig)

rins <- lmer(formula = Persist ~ setting + species + inf_age + parity + context_jt +
	(1|animal_id) + (1|group) + (0+context_jt|animal_id), data = test.data.sig)

# Fit "null models" without random intercepts & slope for individual identity
ri.null <- glmer(formula = Persist ~ setting + species + inf_age + parity + context_jt +
	(1|group), family = binomial, data = test.data.sig, control = contr)
	
ris.null <- glmer(formula = Persist ~ setting + species + inf_age + parity + context_jt +
	(1|group), family = binomial, data = test.data.sig, control = contr)

summary(ri)
AIC(ri,ris)
# df      AIC
# ri   8 611.2300
# ris 11 598.1605

as.data.frame(anova(ri.null, ri, test="Chisq"))
#         npar      AIC      BIC    logLik deviance    Chisq Df  Pr(>Chisq)
# ri.null    7 619.4657 650.8046 -302.7329 605.4657       NA NA          NA
# ri         8 611.2300 647.0458 -297.6150 595.2300 10.23575  1 0.001377447

as.data.frame(anova(ris.null, ris, test="Chisq"))
#          npar      AIC      BIC    logLik deviance    Chisq Df   Pr(>Chisq)
# ris.null    7 619.4657 650.8046 -302.7329 605.4657       NA NA           NA
# ris        11 598.1605 647.4072 -288.0803 576.1605 29.30522  4 6.777033e-06


# Calculate repeatability using rptR ---
#rep_ri <- rpt(Persist ~ setting + species + inf_age + parity + context_jt +
#	(1|animal_id) + (1|group), grname = c("animal_id", "group"), data = test.data.sig, datatype = "Binary", nboot = 1000, 
#	npermut = 0, adjusted = FALSE)
#print(rep_ri)
# Repeatability estimation using the glmm method and logit link 
#
# Repeatability for animal_id
# --------------------------------
#   Link-scale approximation:
#   R  = 0.077
# SE = 0.039
# CI = [0, 0.151]
# P  = 0.000637 [LRT]
# NA [Permutation]

# Original-scale approximation:
#   R  = 0.072
# SE = 0.036
# CI = [0, 0.134]
# P  = 0.000637 [LRT]
# NA [Permutation]

# Repeatability for group
# --------------------------------
#   Link-scale approximation:
#   R  = 0.07
# SE = 0.04
# CI = [0, 0.14]
# P  = 0.274 [LRT]
# NA [Permutation]
# 
# Original-scale approximation:
#   R  = 0.066
# SE = 0.036
# CI = [0, 0.127]
# P  = 0.274 [LRT]
# NA [Permutation]


# Plot individual differences ---
randomSims <- REsim(ri, n.sims = 1000)
head(randomSims[randomSims$groupFctr=="animal_id",])
head(randomSims[randomSims$groupFctr=="group",])

# Add info on the individuals
randomSims <- merge(randomSims[randomSims$groupFctr=="animal_id",],
	test.data.sig[!duplicated(test.data.sig$animal_id),c("animal_id","setting","species","specset")],
	by.x = "groupID",by.y="animal_id")

# Add identifier to color individuals uniquely
randomSims$ID <- ifelse(randomSims$groupID %in% c("CAH","FRI","KER","MAN"), randomSims$groupID, "Other Individuals")
randomSims$COL <- ifelse(randomSims$specset=="Borcaptive",col.born.capt,ifelse(randomSims$specset=="Borwild",col.born.wild,ifelse(randomSims$specset=="Sumcaptive",col.suma.capt,ifelse(randomSims$specset=="Sumwild",col.suma.wild,NA))))
cols <- c("Borcaptive" = col.born.capt, "Borwild" = col.born.wild, "Sumcaptive" = col.suma.capt, "Sumwild" = col.suma.wild)

# Add population intercept and site specific differences for easier interpretation of realized sequence use
randomSims[randomSims$species == "Bor",]$mean <- randomSims[randomSims$species == "Bor",]$mean + fixef(ri)["(Intercept)"]
randomSims[randomSims$species == "Sum",]$mean <- randomSims[randomSims$species == "Sum",]$mean + (fixef(ri)["(Intercept)"] + fixef(ri)["speciesSum"])

# Sort data by value and specset
table(randomSims$specset)
sort.order <- c("Borcaptive","Borwild","Sumcaptive","Sumwild")
randomSims <- arrange(randomSims, mean)
randomSims <- arrange(randomSims, factor(specset, levels = sort.order))

# Tell ggplot that you have an ordered factor already
randomSims$groupID <- factor(randomSims$groupID, levels = randomSims$groupID) #for ordering: [order(randomSims$mean)]

#svg(filename=paste(path,outfile,".svg",sep=""), height=50/25.4, width=70/25.4, family="Arial", pointsize = 9)
ggplot() +
	geom_errorbar(data = randomSims, aes(x = groupID, ymin = plogis(mean-sd), ymax = plogis(mean+sd), col=specset)) +
	geom_point(data = randomSims, aes(x = groupID, y = plogis(mean), col=specset), size = 2) +
#	geom_vline(xintercept = c(4.5,12.5,17.5), linetype=2) + 
	theme_classic() +
	ylab("BT gestural redoings") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#	annotate("text",x = 3, y = 0.4, label = "R = 0.1", size = 4) +
	scale_color_manual(values = cols, name="Species-setting class", labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild"))


# -----------------------------------------------------
# Plot output of linear mixed-effects models (lmer) ---
test.data.sig$context_jt <- as.factor(test.data.sig$context_jt)
test.data.sig$context_jt <- relevel(test.data.sig$context_jt, "1")

# Plot model predictions and raw data using lmer models (individual variation in plasticity)
# augment() throws an error with a glmm. However, we only use the ".fitted" column, which is not affected by that warning.
# I tested this by using x <- predict(ri) ; cor(x,RI$.fitted) ; plot(x,RI$.fitted) # --> values are the same
RI <- augment(ri) %>%
	dplyr::select(Persist, context_jt , animal_id, .fitted)
# augment() throws an error with a glmm. However, we only use the ".fitted" column, which is not affected by that warning.
# I tested this by using x <- predict(ris) ; cor(x,RIS$.fitted) ; plot(x,RIS$.fitted) # --> values are the same
RIS <- augment(ris) %>%
	dplyr::select(Persist, context_jt , animal_id, .fitted)
RI$.fittedRIS <- RIS$.fitted
df <- RI %>%
	dplyr::group_by(animal_id, context_jt) %>%
	dplyr::summarise(RI = mean(.fitted), RIS = mean(.fittedRIS))%>%
	gather(type, Value, `RI`:`RIS`)
df$animal_id <- as.character(df$animal_id)
df <- merge(df, test.data.sig[!duplicated(test.data.sig$animal_id),c("animal_id","specset")], by.x = "animal_id",by.y="animal_id")
df$ID <- ifelse(df$animal_id %in% c("CAH","FRI","KER","MAN"), df$animal_id, "Other Individuals")
df$COL <- ifelse(df$specset=="Borcaptive",col.born.capt,ifelse(df$specset=="Borwild",col.born.wild,ifelse(df$specset=="Sumcaptive",col.suma.capt,ifelse(df$specset=="Sumwild",col.suma.wild,NA))))

plot_ri <- ggplot(df[df$type == "RI",], aes(x = context_jt, y = plogis(Value), group = animal_id, color = specset)) +
	geom_line() +
	theme_classic() +
	labs(y="Gestural redoings", x="") +
	ggtitle("Random Intercept") +
	ylim(0,1) +
	scale_color_manual(values = cols, name="Species-setting class", labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild")) +
	scale_x_discrete("Social context", limits = c("1","0"), labels = c("co-locomotion","other")) +
	guides(color = guide_legend(nrow = 1, byrow = TRUE))
	
plot_ris <- ggplot(df[df$type == "RIS",], aes(x = context_jt, y = plogis(Value), group = animal_id, color = specset)) +
	geom_line() +
	theme_classic() +
	labs(y="", x="") +
	ggtitle("Random Intercept\n and Slope") +
	ylim(0,1) +
	scale_x_discrete("Social context", limits = c("1","0"), labels = c("co-locomotion","other")) +
	scale_color_manual(values = cols, labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild"))
mylegend <- g_legend(plot_ri)
grid.arrange(arrangeGrob(plot_ri + theme(legend.position="none"),
	plot_ris + theme(legend.position="none"), ncol=2), mylegend, heights=c(10, 2))

# The random regression looks a bit weird because connecting lines are not parallel. This is because we back-transformed the data.
# On the logit-scale, connecting lines are parallel (i.e. only the intercepts vary).




# --------------------------------------------------------------------------------------------------------
# (6) Analyses: 2nd variable: Responsiveness to infant requests ------------------------------------------
# ----------------------------------------------------
# (6.1) Frequentist approach -------------------------
# Random intercept model
ri_2 <- glmer(formula = ASO ~ setting + species + inf_age + parity + context_fs +
	(1|ID_rec) + (1|group), family = binomial, data = test.data.out, control = contr)

# Random slope model
ris_2 <- glmer(formula = ASO ~ setting + species + inf_age + parity + context_fs +
	(1|ID_rec) + (0+context_fs|ID_rec) + (1|group), family = binomial, data = test.data.out, control = contr)

# Check models with lmer
rin_2 <- lmer(formula = ASO ~ setting + species + inf_age + parity + context_fs +
	(1|ID_rec) + (1|group), data = test.data.out)

rins_2 <- lmer(formula = ASO ~ setting + species + inf_age + parity + context_fs +
	(1|ID_rec) + (0+context_fs|ID_rec) + (1|group), data = test.data.out)

length(residuals(ri_2)) #3446
length(residuals(ris_2)) #3446

# Fit "null models" without random intercepts & slope for individual identity
ri_2.null <- glmer(formula = ASO ~ setting + species + inf_age + parity + context_fs +
	(1|group), family = binomial, data = test.data.out, control = contr)
ris_2.null <- glmer(formula = ASO ~ setting + species + inf_age + parity + context_fs +
	(1|group), family = binomial, data = test.data.out, control = contr)

summary(ri_2)
AIC(ri_2,ris_2)
# ri_2   8 3915.324
# ris_2  9 3828.484

as.data.frame(anova(ris_2.null, ris_2, test="Chisq"))
#            npar      AIC      BIC    logLik deviance    Chisq Df   Pr(>Chisq)
# ris_2.null    7 3928.389 3971.404 -1957.194 3914.389       NA NA           NA
# ris_2         9 3828.484 3883.789 -1905.242 3810.484 103.9044  2 2.738106e-23

as.data.frame(anova(ri_2.null, ri_2, test="Chisq"))
#           npar      AIC      BIC    logLik deviance    Chisq Df   Pr(>Chisq)
# ri_2.null    7 3928.389 3971.404 -1957.194 3914.389       NA NA           NA
# ri_2         8 3915.324 3964.484 -1949.662 3899.324 15.06486  1 0.0001038793

#Calculate repeatability using rptR
rep_ri2 <- rpt(ASO ~ setting + species + inf_age + parity + context_fs +
	(1|animal_id) + (1|group), grname = c("animal_id","group"), data = test.data.out, datatype = "Binary", nboot = 1000,
	npermut = 0, adjusted = FALSE)
print(rep_ri2)
# Repeatability estimation using the glmm method and logit link 
#
# Repeatability for animal_id
# --------------------------------
#   Link-scale approximation:
#   R  = 0.017
# SE = 0.009
# CI = [0, 0.036]
# P  = 4.42e-05 [LRT]
# NA [Permutation]
# 
# Original-scale approximation:
#   R  = 0.015
# SE = 0.009
# CI = [0, 0.034]
# P  = 4.42e-05 [LRT]
# NA [Permutation]
# 
# Repeatability for group
# --------------------------------
#   Link-scale approximation:
#   R  = 0.114
# SE = 0.049
# CI = [0, 0.173]
# P  = 1.02e-05 [LRT]
# NA [Permutation]
# 
# Original-scale approximation:
#   R  = 0.102
# SE = 0.045
# CI = [0, 0.159]
# P  = 1.02e-05 [LRT]
# NA [Permutation]


# Plot individual differences ---
randomSims <- REsim(ri_2, n.sims = 1000)
head(randomSims[randomSims$groupFctr=="ID_rec",])
head(randomSims[randomSims$groupFctr=="group",])

# Add info on the individuals
randomSims <- merge(randomSims[randomSims$groupFctr=="ID_rec",],
                    test.data.out[!duplicated(test.data.out$ID_rec),c("ID_rec","setting","species","specset")],
                    by.x = "groupID",by.y="ID_rec")

# Add identifier to color individuals uniquely
randomSims$ID <- ifelse(randomSims$groupID %in% c("CAH","FRI","KER","MAN"), randomSims$groupID, "Other Individuals")
randomSims$COL <- ifelse(randomSims$specset=="Borcaptive",col.born.capt,ifelse(randomSims$specset=="Borwild",col.born.wild,ifelse(randomSims$specset=="Sumcaptive",col.suma.capt,ifelse(randomSims$specset=="Sumwild",col.suma.wild,NA))))
cols <- c("Borcaptive" = col.born.capt, "Borwild" = col.born.wild, "Sumcaptive" = col.suma.capt, "Sumwild" = col.suma.wild)

# Add population intercept and site specific differences for easier interpretation of realized sequence use
randomSims[randomSims$species == "Bor",]$mean <- randomSims[randomSims$species == "Bor",]$mean + fixef(ri_2)["(Intercept)"]
randomSims[randomSims$species == "Sum",]$mean <- randomSims[randomSims$species == "Sum",]$mean + (fixef(ri_2)["(Intercept)"] + fixef(ri_2)["speciesSum"])

# Sort data by value and specset
table(randomSims$specset)
sort.order <- c("Borcaptive","Borwild","Sumcaptive","Sumwild")
randomSims <- arrange(randomSims, mean)
randomSims <- arrange(randomSims, factor(specset, levels = sort.order))

# Tell ggplot that you have an ordered factor already
randomSims$groupID <- factor(randomSims$groupID, levels = randomSims$groupID) #for ordering: [order(randomSims$mean)]

#svg(filename=paste(path,outfile,".svg",sep=""), height=50/25.4, width=70/25.4, family="Arial", pointsize = 9)
ggplot() +
  geom_errorbar(data = randomSims, aes(x = groupID, ymin = plogis(mean-sd), ymax = plogis(mean+sd), col=specset)) +
  geom_point(data = randomSims, aes(x = groupID, y = plogis(mean), col=specset), size = 2) +
  #	geom_vline(xintercept = c(4.5,12.5,17.5), linetype=2) + 
  theme_classic() +
  ylab("BT responsiveness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  #	annotate("text",x = 3, y = 0.4, label = "R = 0.1", size = 4) +
  scale_color_manual(values = cols, name="Species-setting class", labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild"))


# Plot model predictions and raw data using lmer models (individual variation in plasticity)
# augment() throws an error with a glmm. However, we only use the ".fitted" column, which is not affected by that warning.
# I tested this by using x <- predict(ri) ; cor(x,RI$.fitted) ; plot(x,RI$.fitted) # --> values are the same
RI <- augment(ri_2) %>%
  dplyr::select(ASO, context_fs , ID_rec, .fitted)
# augment() throws an error with a glmm. However, we only use the ".fitted" column, which is not affected by that warning.
# I tested this by using x <- predict(ris) ; cor(x,RIS$.fitted) ; plot(x,RIS$.fitted) # --> values are the same
RIS <- augment(ris_2) %>%
  dplyr::select(ASO, context_fs , ID_rec, .fitted)
RI$.fittedRIS <- RIS$.fitted
df <- RI %>%
  dplyr::group_by(ID_rec, context_fs) %>%
  dplyr::summarise(RI = mean(.fitted), RIS = mean(.fittedRIS))%>%
  gather(type, Value, `RI`:`RIS`)
df$ID_rec <- as.character(df$ID_rec)
df <- merge(df, test.data.out[!duplicated(test.data.out$ID_rec),c("ID_rec","specset")], by.x = "ID_rec",by.y="ID_rec")
df$ID <- ifelse(df$ID_rec %in% c("CAH","FRI","KER","MAN"), df$ID_rec, "Other Individuals")
df$COL <- ifelse(df$specset=="Borcaptive",col.born.capt,ifelse(df$specset=="Borwild",col.born.wild,ifelse(df$specset=="Sumcaptive",col.suma.capt,ifelse(df$specset=="Sumwild",col.suma.wild,NA))))

plot_ri <- ggplot(df[df$type == "RI",], aes(x = context_fs, y = plogis(Value), group = ID_rec, color = specset)) +
  geom_line() +
  theme_classic() +
  labs(y="Responsiveness", x="") +
  ggtitle("Random Intercept") +
  ylim(0,1) +
  scale_color_manual(values = cols, name="Species-setting class", labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild")) +
  scale_x_discrete("Social context", limits = c("1","0"), labels = c("begging","other")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

plot_ris <- ggplot(df[df$type == "RIS",], aes(x = context_fs, y = plogis(Value), group = ID_rec, color = specset)) +
  geom_line() +
  theme_classic() +
  labs(y="", x="") +
  ggtitle("Random Intercept\n and Slope") +
  ylim(0,1) +
  scale_x_discrete("Social context", limits = c("1","0"), labels = c("begging","other")) +
  scale_color_manual(values = cols, name="Species-setting class", labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild"))
mylegend <- g_legend(plot_ri)
grid.arrange(arrangeGrob(plot_ri + theme(legend.position="none"),
                         plot_ris + theme(legend.position="none"), ncol=2), mylegend, heights=c(10, 2))


# --------------------------------------------------------------------------------------------------------
# (7) Save workspace -------------------------------------------------------------------------------------
save.image(paste0(path,"orang_motoff_MF.RData"))