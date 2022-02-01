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
path <- "G:\\Knief\\Documents\\CoAuthoredManuscripts\\OrangUtan_RepertoireSimilarityPlasticity\\"
dat <- read.table(paste0(path,"Data\\orangutanR_inteng_repair.csv"), header=TRUE, sep=",", stringsAsFactors=TRUE)

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

as.data.frame(anova(ris.null, ris, test="Chisq")) 											# -> I PUT IN THE NEW VALUES HERE
#          npar      AIC      BIC    logLik deviance    Chisq Df   Pr(>Chisq)
# ris.null    7 619.4657 650.8046 -302.7329 605.4657       NA NA           NA
# ris        11 598.1605 647.4072 -288.0803 576.1605 29.30522  4 6.777033e-06

# Plot population-level effect of context on sequence use ---
ggplot(test.data.sig,
	aes(x = context_jt, y = Persist, color = species)) +
	geom_jitter(size=0.5) +
	geom_point(aes(y = plogis(predict(ri,re.form = NA))),size=1) +
	theme_classic() +
	ylab("Persistence")#+
#	scale_x_discrete("Context", limits = c("non-travel", "travel"))

# Plot individual differences in persistence in travel versus non-travel context ---
pr <- data.frame(animal_id = test.data.sig$animal_id, context = test.data.sig$context_jt)
pr$fit <- plogis(predict(ri))

ggplot() +
	geom_line(data = pr, aes(x = context, y = pr[,"fit"], group = animal_id, color = animal_id), size=0.3) +
	labs(y = "Predicted persistence",color = "animal_id") +
	theme_classic() +
	theme(legend.position="none") +
	guides(color = guide_legend(nrow = 4))#+
#	scale_x_discrete("Context", limits = c("non-travel context", "travel"))


# Calculate repeatability ---
print(VarCorr(ri),comp = c("Variance","Std.Dev."))
# Groups    Name        Variance Std.Dev.
# animal_id (Intercept) 0.57792  0.76021 
# group     (Intercept) 0.54589  0.73884 

VarCorr(ri)$"animal_id"[1] / (VarCorr(ri)$"animal_id"[1] + VarCorr(ri)$"group"[1] + attr(VarCorr(ri), "sc")^2)
# [1] 0.2721161

print(VarCorr(rin),comp = c("Variance","Std.Dev."))
# Groups    Name        Variance Std.Dev.
# animal_id (Intercept) 0.016028 0.12660 
# group     (Intercept) 0.011995 0.10952 
# Residual              0.141335 0.37595  

VarCorr(rin)$"animal_id"[1] / (VarCorr(rin)$"animal_id"[1] + VarCorr(rin)$"group"[1] + attr(VarCorr(rin), "sc")^2)
# [1] 0.09464031

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

# Simulate model (lmer) to get a posterior distribution for all variance components ---
set.seed(1)
simulated <- sim(rin, n.sim = 1000)
posterior_animal_id <- apply(simulated@ranef$"animal_id"[ , , 1],1,var)
posterior_group <- apply(simulated@ranef$"group"[ , , 1],1,var)
posterior_residual <- simulated@sigma^2

quantile(posterior_animal_id / (posterior_animal_id + posterior_residual), prob=c(0.025, 0.5, 0.975))
#       2.5%        50%      97.5% 
# 0.05904336 0.10159168 0.16208661  

# Coefficients of variance: dividing the square root of the among-individual variance by the intercept (i.e. mean trait value) ---
CVi <- sqrt(posterior_animal_id) / summary(rin)$coefficients[1]
quantile(CVi,prob=c(0.025, 0.5, 0.975))
# 2.5%      50%       97.5% 
# 0.7998096 1.0966744 1.4263180 

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
	scale_color_manual(values = cols, name="Species-setting-class", labels=c("Bornean captive","Bornean wild","Sumatran captive","Sumatran wild"))
mylegend <- g_legend(plot_ri)
grid.arrange(arrangeGrob(plot_ri + theme(legend.position="none"),
	plot_ris + theme(legend.position="none"), ncol=2), mylegend, heights=c(10, 2))

# The random regression looks a bit weird because connecting lines are not parallel. This is because we back-transformed the data.
# On the logit-scale, connecting lines are parallel (i.e. only the intercepts vary).


# ----------------------------------------------------
# (5.2) Bayesian approach ----------------------------
# WORKS ONLY IF FAMILY TERM IS NOT INCLUDED, AND THUS IGNORING THAT THE DATA IS BINARY
my.cores <- detectCores()
m1_brm <- brm(Persist ~ setting + species + inf_age + parity + context_jt +
	(1|animal_id) + (1|group),
	data = test.data.sig,
#	family = bernoulli(link = "logit"),
	warmup = 500,
	iter = 3000,
	thin=2,
	chains = 2,
	inits = "random",
	cores = my.cores,
	seed = 12345)
m1_brm <- add_criterion(m1_brm, "waic")

m1_brm <- readRDS("m1_brm.rds")
summary(m1_brm)

colnames(posterior_samples(m1_brm))[1:8]

# Calculate repeatability ---
var.animal_id <- posterior_samples(m1_brm)$"sd_animal_id__Intercept"^2
var.group <- posterior_samples(m1_brm)$"sd_group__Intercept"^2
var.res <- posterior_samples(m1_brm)$"sigma"^2

RPer <- var.animal_id / (var.animal_id + var.group + var.res)
mean(RPer) ; HPDinterval(as.mcmc(RPer),0.95)
# [1] 0.1085193
# lower     upper
# var1 0.01010544 0.2285061
# attr(,"Probability")
# [1] 0.95

CVi <- sqrt(var.animal_id) / mean(test.data.sig$Persist)
mean(CVi) ; HPDinterval(as.mcmc(CVi),0.95)
# [1] 0.6920443
# lower    upper
# var1 0.3456064 1.118442
# attr(,"Probability")
# [1] 0.95

# Plot brms output ---
# I HAVE TROUBLES TO REENACT THE FIRST SECTION OF THIS CODE - WHAT DO THEY "SEPARATE" HERE?
posteriorBT <- posterior_samples(m1_brm)[,10:32] %>%
	gather(animal_id, value, "r_animal_id[CAH,Intercept]" : "r_animal_id[WAT,Intercept]")%>%
	separate(animal_id, c(NA,NA,NA,"animal_id",NA), sep = "([\\_\\[\\,])", fill = "right") %>%
	left_join(select(data[!duplicated(data$animal_id),], animal_id, species))
posteriorBT[posteriorBT$species == "Bor",]$value <- posteriorBT[posteriorBT$species == "Bor",]$value + fixef(m1_brm, pars = "Intercept")[1]
posteriorBT[posteriorBT$species == "Sum",]$value <- posteriorBT[posteriorBT$species == "Sum",]$value + fixef(m1_brm, pars = "Intercept")[1] + fixef(m1_brm, pars = "speciesSum")[1]
posteriorBT$col <- ifelse(posteriorBT$animal_id %in% c("CAH","KER","FRI","MAN"), posteriorBT$animal_id, "Other individuals")
posteriorBT <- posteriorBT %>%
	dplyr::group_by(animal_id) %>%
	dplyr::mutate(meanBT = mean(value))%>%
	dplyr::ungroup()

BT <- ggplot()+
	ggridges::geom_density_ridges(data = posteriorBT, aes(x = value, y = reorder(as.factor(animal_id), meanBT), height = ..density.., fill = col,scale = 3), alpha = 0.6) +
	geom_point(data = posteriorBT[!duplicated(posteriorBT$animal_id),], aes(x = meanBT, y = as.factor(animal_id), col = species), size = 1) +
	labs(y = "", x = "BT persistence", fill = "ID") +
	theme_classic() +
	scale_fill_manual(values = c("#F8766D","#C77CFF","#7CAE00","#FFCC00","gray"))

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
	(1|group), family = binomial, data = test.data.out, control = contr)					# -> ri_2.null AND ris_2.null ARE IDENTICAL

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

# Plot population-level effect of context on responsiveness ---
ggplot(test.data.out,
	aes(x = context_fs, y = ASO, color = species)) +
	geom_jitter(size=0.5) +
	geom_point(aes(y = plogis(predict(ri_2,re.form = NA))),size=1) +
	theme_classic() +
	ylab("Responsiveness")#+
#	scale_x_discrete("Context", limits = c("not food-sharing", "food-sharing"))

# Plot individual differences in responsiveness in food-sharing versus non-food-sharing interactions ---
pr <- data.frame(animal_id = test.data.out$ID_rec, context = test.data.out$context_fs)
pr$fit <- plogis(predict(ri_2))
ggplot() +
	geom_line(data = pr, aes(x = context, y = pr[,"fit"], group = animal_id, color = animal_id), size=0.3) +
	labs(y = "Predicted responsiveness",color = "animal_id") +
	theme_classic() +
	theme(legend.position="none") +
	guides(color = guide_legend(nrow = 4))#+
#	scale_x_discrete("Context", limits = c("not food sharing", "food sharing"))

# Calculate repeatability ---
print(VarCorr(ri_2),comp = c("Variance","Std.Dev."))
# Groups Name        Variance Std.Dev.
# ID_rec (Intercept) 0.08267  0.28752 
# group  (Intercept) 0.56671  0.75280 

VarCorr(ri_2)$"ID_rec"[1] / (VarCorr(ri_2)$"group"[1] + attr(VarCorr(ri_2), "sc")^2) 			# -> I THINK THAT YOU NEED TO ADD VarCorr(ri_2)$"ID_rec"[1] ALSO TO THE DENOMINATOR!!!
# [1] 0.0532948

# Calculate repeatability lmer ---
print(VarCorr(rin_2),comp = c("Variance","Std.Dev."))
# Groups   Name        Variance  Std.Dev.
# ID_rec   (Intercept) 0.0038991 0.062443
# group    (Intercept) 0.0428306 0.206956
# Residual             0.1894890 0.435303  

VarCorr(rin_2)$"animal_id"[1] / (VarCorr(rin_2)$"group"[1] + attr(VarCorr(rin_2), "sc")^2) 		# -> I THINK THAT YOU NEED TO USE "ID_rec" AND THAT YOU NEED TO ADD VarCorr(ri_2)$"ID_rec"[1] ALSO TO THE DENOMINATOR!!!

#Calculate repeatability using rptR
rep_ri2 <- rpt(ASO ~ setting + species + inf_age + parity + context_fs +
	(1|animal_id) + (1|group), grname = c("animal_id","group"), data = test.data.out, datatype = "Binary", nboot = 1000, 	# -> I THINK THAT YOU NEED TO USE "ID_rec"
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

# Simulate model (lmer) to get a posterior distribution for all variance components ---
set.seed(1)
simulated <- sim(rin_2, n.sim = 1000)
posterior_animal_id <- apply(simulated@ranef$"ID_rec"[ , , 1],1,var)
posterior_group <- apply(simulated@ranef$"group"[ , , 1],1,var)
posterior_residual <- simulated@sigma^2

quantile(posterior_animal_id / (posterior_animal_id + posterior_residual), prob=c(0.025, 0.5, 0.975))
# 2.5%        50%      97.5% 
# 0.01095754 0.01905934 0.03007454

# Coefficients of variance: dividing the square root of the among-individual variance by the intercept (i.e. mean trait value) ---
CVi <- sqrt(posterior_animal_id) / summary(rin_2)$coefficients[1]
quantile(CVi,prob=c(0.025, 0.5, 0.975))
# 2.5%       50%     97.5% 
# 0.1044580 0.1364462 0.1713143   

# Plot individual differences ---
randomSims <- REsim(ri_2, n.sims = 1000)
head(randomSims[randomSims$groupFctr=="ID_rec",])

# Add the species of the individual
randomSims <- merge(randomSims[randomSims$groupFctr=="ID_rec",],
	test.data.out[!duplicated(test.data.out$ID_rec),c("ID_rec","species")],
	by.x = "groupID",by.y="ID_rec")
# Add identifier to color individuals uniquely
randomSims$motherID <- ifelse(randomSims$groupID %in% c("CAH","KER","FRI","MAN"), randomSims$groupID, "Other")

# Add population intercept and site specific differences for easier interpretation of realized sequence use
randomSims[randomSims$species == "Bor",]$mean <- randomSims[randomSims$species == "Bor",]$mean + fixef(ri_2)["(Intercept)"]
randomSims[randomSims$species == "Sum",]$mean <- randomSims[randomSims$species == "Sum",]$mean + (fixef(ri_2)["(Intercept)"] + fixef(ri_2)["speciesSum"])

# Sort data
randomSims$groupID <- factor(randomSims$groupID, levels = randomSims$groupID) #for odering: [order(randomSims$mean)]
ggplot()+
	geom_errorbar(data = randomSims, aes(x = groupID, ymin = plogis(mean-sd), ymax = plogis(mean+sd), color = species)) +
	geom_point(data = randomSims, aes(x = groupID, y = plogis(mean), fill = motherID), shape = 21, size = 2) +
	theme_classic() +
	ylab("BT responsiveness") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
#	annotate("text",x = 3, y = 0.4, label = "R = 0.05", size = 4) +
	scale_fill_manual(values = c("#F8766D","#C77CFF","#7CAE00","#FFCC00","gray"))

# -----------------------------------------------------
# Plot output of linear mixed-effects models (lmer) ---
test.data.out$context_fs <- as.factor(test.data.out$context_fs)
test.data.out$context_fs <- relevel(test.data.out$context_fs, "1")
test.data.out$context_fs <- as.numeric(test.data.out$context_fs)

# Plot model predictions and raw data using lmer models (individual variation in plasticity) ---
RI<-augment(rin_2) %>%
	dplyr::select(ASO, context_fs , ID_rec, .fitted)
RIS<-augment(rins_2) %>%
	dplyr::select(ASO, context_fs , ID_rec, .fitted)
RI$.fittedRIS<-RIS$.fitted
df <- RI %>%
	dplyr::group_by(ID_rec, context_fs) %>%
	dplyr::summarise(RI = mean(.fitted), RIS = mean(.fittedRIS))%>%
	gather(type, Value, `RI`:`RIS`)
df$ID_rec <- as.character(df$ID_rec)
df$ID <- ifelse(df$ID_rec %in% c("CAH","KER","FRI","MAN"), df$ID_rec, "Other Individuals")

plot_ri <- ggplot(df[df$type == "RI",], aes(x = context_fs, y = Value, group = ID_rec,color = ID)) +
	geom_line() +
	theme_classic() +
	labs(y="Responsiveness", x="") +
	ggtitle("(a) Random Intercept") +
	ylim(0,1) +
#	scale_y_continuous("", breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
	scale_color_manual(values = c("#F8766D","#C77CFF","#7CAE00","#00BFC4","gray")) +
	scale_x_discrete("Social context", limits = c("1", "0"), labels = c("begging","other")) +
	guides(color = guide_legend(nrow = 1,byrow = TRUE))

plot_ris <- ggplot(df[df$type == "RIS",], aes(x = context_fs, y = Value, group = ID_rec,color = ID)) +
	geom_line() +
	theme_classic() +
	labs(y="", x="") +
	ggtitle("(b) Random Intercept and Slope") +
	ylim(0,1) +
#	scale_y_continuous("", breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
	scale_x_discrete("Social context", limits = c("1","0"), labels = c("begging","other")) +
	scale_color_manual(values = c("#F8766D","#C77CFF","#7CAE00","#00BFC4","gray"))
	mylegend <- g_legend(plot_ri)
grid.arrange(arrangeGrob(plot_ri + theme(legend.position="none"),
	plot_ris + theme(legend.position="none"), ncol=2), mylegend, heights=c(10, 2))

# --------------------------------------------------------------------------------------------------------
# (7) Save workspace -------------------------------------------------------------------------------------
save.image(paste0(path,"orang_motoff_MF.RData"))