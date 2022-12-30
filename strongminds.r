library(rstan)
library(cmdstanr)
library(tidyverse)

setwd("/users/mattlerner/desktop/evaluations/strongminds")


mc.cores = parallel::detectCores()

################### Read in StrongMinds programmes data ###################

strongminds.data.1 <- read.csv("data/hli/StrongMinds data via HLI - Sheet1.csv")
colnames(strongminds.data.1) <- c('program','type','country','year','pre_n','pre_mean','pre_sd','post_n','post_mean','post_sd','mean_change','change_relative_core')
strongminds.data.1$pre_n <- as.numeric(gsub(",","",strongminds.data.1$pre_n))
strongminds.data.1$post_n <- as.numeric(gsub(",","",strongminds.data.1$post_n))
strongminds.data.1$attrition <- strongminds.data.1$post_n - strongminds.data.1$pre_n
strongminds.data.1$pct_attrition <- strongminds.data.1$attrition/strongminds.data.1$pre_n


################### Read in HLI's estimates for StrongMinds program effects ###################

hli.strongminds.effects.raw <- read.csv("data/hli/SM's Point Estimate CEA - Cleaner Point Estimate CEA.csv")
hli.strongminds.effects <- hli.strongminds.effects.raw[1:8,1:14]

################### Clean up frames ###################
hli.strongminds.effects$variable.cost <- as.numeric(gsub("\\$","",hli.strongminds.effects$Avg..Variable.Cost.Per.Person))
hli.strongminds.effects$program.budget.pct <- as.numeric(gsub("\\%","",hli.strongminds.effects$Percent.of.Program.Budget))/100


################### Calculate the program pre-participation mean PHQ-9 score for different countries and programs ###################

# weighted mean
sm.program.level <- strongminds.data.1 %>%
  group_by(country, type) %>%
  summarise(pre.mean = weighted.mean(pre_mean, pre_n),
            pre.sd = sqrt(weighted.mean(pre_sd**2, pre_n))
            )

################### Functions and parameters ###################
N <- 100000

# inverse quantile function for exponential distribution
invq.exp <- function(p, quantile) {
  lambda <- (-1 * (log(1-p)))/quantile
  return(lambda)
}


cohen.d <- function(mu1, mu2, se1, se2, n1, n2) {
  # back out sds from ses
  s1 <- se1*sqrt(n1)
  s2 <- se2*sqrt(n2)
  
  s.numerator <- ((n1 - 1)*(s1^2)) + ((n2 - 1)*(s2^2))
  s.denominator <- n1 + n2 - 2
  s <- sqrt(s.numerator/s.denominator)
  output <- (mu1 - mu2)/s
  return(output)
}

# you can derive this just using the definition of mean and standard deviation
# the idea here is to imagine that, for all of the attriters, the treatment
# made approximately no difference, e.g. that all attriters had the mean treatment
# value before and after the intervention.
deflate.mean.and.sd <- function(pre_mean, pre_sd, post_mean, post_sd, pre_n, post_n) {
  attrition <- pre_n - post_n
  deflated.mean <- ((post_mean*post_n) + (attrition*pre_mean))/pre_n
  deflated.sd <- sqrt(1/pre_n) * sqrt((post_n * (post_sd^2)) + attrition*(pre_mean-deflated.mean)^2)
  return(list(mean = deflated.mean, sd = deflated.sd))
}

################### Background work ###################

######## Back out distributions for country-level PHQ-9 distributions on the assumption they are roughly exponential
# See e.g. this reference https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6115508/
# See also PHQ-9 instruction manual https://www.pcpcc.org/sites/default/files/resources/instructions.pdf

# Share with depression by country via Our World in Data
# https://ourworldindata.org/grapher/share-with-depression?tab=table
# Per this article https://pubmed.ncbi.nlm.nih.gov/31439359/
# the vast majority of depression cases are major depression, not dysthymia

# Per this study: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3281183/
# We use a cutoff score of 10 for major depression (recommendation of 8-11)

depression.cutoff <- 10

### Uganda, 4.92%
# So we want parameters for an exponential distribution such that 4.92% of mass lies above 10
# e.g. the 95.08th percentile is 10, the depression cutoff
uganda.parameter <- invq.exp(0.9508, depression.cutoff)
uganda.phq.distribution <- rexp(N, rate=uganda.parameter)
uganda.mean <- mean(uganda.phq.distribution)

### Zambia, 3.64%
zambia.parameter <- invq.exp(0.9636, depression.cutoff)
zambia.phq.distribution <- rexp(N, rate=zambia.parameter)
zambia.mean <- mean(zambia.phq.distribution)

# distributions for depressed people
uganda.depressed.distribution <- uganda.phq.distribution[uganda.phq.distribution >= depression.cutoff]
zambia.depressed.distribution <- zambia.phq.distribution[zambia.phq.distribution >= depression.cutoff]

distros <- list(Zambia = zambia.depressed.distribution, Uganda = uganda.depressed.distribution)

################### ################### ################### ################### ################### 
# Part 1: Pressure-testing our first evaluation to see if linearizing disability weights made a big difference.
################### ################### ################### ################### ################### 

################### Calculating change in incidence of depression for a given program ###################

# Fixing merge
hli.strongminds.effects[hli.strongminds.effects$Type == "Covid, Refugee & Partner (Peer)", "Type"] <- "Refugee-Hybrid"
programs.x.effects <- merge(sm.program.level, hli.strongminds.effects, by.x=c("country","type"), by.y=c("Country","Type"))

# we assume an exponential distribution of scores -- the mean and the SD are the same
programs.x.effects$effect <- as.numeric(programs.x.effects$Total.Effect) * programs.x.effects$pre.sd

# now we want to know what the distribution of participants in each depression category was
# PRIOR to the intervention

# See PHQ-9 instruction manual https://www.pcpcc.org/sites/default/files/resources/instructions.pdf
mild.depression.cutoff <- 5
moderate.depression.cutoff <- 10
severe.depression.cutoff <- 15 # manual has cutoffs at 15 for "moderately severe" and 20 for "severe". We combine them here.

# percentage-point incidence in the incidence of mild depression
get.mild.depression.incidence <- function(country, offset) {
  phq9.dist <- distros[[country]] - offset
  normal.mild.incidence <- mean(phq9.dist >= mild.depression.cutoff & phq9.dist < moderate.depression.cutoff)
  return(normal.mild.incidence)
}

# percentage-point incidence in the incidence of moderate depression
get.moderate.depression.incidence <- function(country, offset) {
  phq9.dist <- distros[[country]] - offset
  normal.moderate.incidence <- mean(phq9.dist >= moderate.depression.cutoff & phq9.dist < severe.depression.cutoff)
  return(normal.moderate.incidence)
}

# percentage-point incidence in the incidence of severe depression
get.severe.depression.incidence <- function(country, offset) {
  phq9.dist <- distros[[country]] - offset
  normal.severe.incidence <- mean(phq9.dist >= severe.depression.cutoff)
  return(normal.severe.incidence)
}

# calculate portions of each depression level before intervention
programs.x.effects$pre.mild.incidence <- mapply(get.mild.depression.incidence, programs.x.effects$country, 0)
programs.x.effects$pre.moderate.incidence <- mapply(get.moderate.depression.incidence, programs.x.effects$country, 0)
programs.x.effects$pre.severe.incidence <- mapply(get.severe.depression.incidence, programs.x.effects$country, 0)

# calculate portions of each depression level AFTER intervention
programs.x.effects$post.mild.incidence <- mapply(get.mild.depression.incidence, programs.x.effects$country, programs.x.effects$effect)
programs.x.effects$post.moderate.incidence <- mapply(get.moderate.depression.incidence, programs.x.effects$country, programs.x.effects$effect)
programs.x.effects$post.severe.incidence <- mapply(get.severe.depression.incidence, programs.x.effects$country, programs.x.effects$effect)

# calculate the pre/post depression levels for each program
programs.x.effects$mild.effect <- programs.x.effects$post.mild.incidence - programs.x.effects$pre.mild.incidence
programs.x.effects$moderate.effect <- programs.x.effects$post.moderate.incidence - programs.x.effects$pre.moderate.incidence
programs.x.effects$severe.effect <- programs.x.effects$post.severe.incidence - programs.x.effects$pre.severe.incidence

################### Calculate expected DALYs per individual treated ###################

# 2019 disability weights for everything
# https://ghdx.healthdata.org/record/ihme-data/gbd-2019-disability-weights

# Mild depression
mild.depression.dalys <- 0.145

# Moderate depression
moderate.depression.dalys <- 0.396

# Severe depression
severe.depression.dalys <- 0.658

programs.x.effects$daly.effect <- (mild.depression.dalys*programs.x.effects$mild.effect) + (moderate.depression.dalys*programs.x.effects$moderate.effect) + (severe.depression.dalys*programs.x.effects$severe.effect)
programs.x.effects$expected.dalys.preserved <- -1 * programs.x.effects$daly.effect
programs.x.effects$program.ce <- programs.x.effects$variable.cost / programs.x.effects$expected.dalys.preserved

################### Calculate expected cost-effectiveness ###################

overall.ce <- sum((programs.x.effects$variable.cost / programs.x.effects$expected.dalys.preserved) * programs.x.effects$program.budget.pct)
overall.ce

# We currently value a DALY at about $630
rating <- 630/overall.ce
rating

# So this method suggests a rating of about 2.5, still clearing the bar, but with some hand-waving in the analysis here.

################### ################### ################### ################### ################### 
# Part 2: Recovering some parts of HLI's analysis
################### ################### ################### ################### ################### 

################### Bolton 2003, replicating HLI calculation of Cohen's D ###################

#bolton.cohen.d <- cohen.d(17.47, 3.55, 1.1, 1.1, 107, 117) # whole group
bolton.cohen.d <- cohen.d(11.59, 2.38, 0.8, 0.75, 107, 117) # all eligible persons - this seems to be the HLI number
bolton.cohen.d

# calculate effect size for CONTROL group, e.g. standardized pre-post difference among the untreated
bolton.control.cohen <- cohen.d(23.65, 21.14, 6.3, 8.19, 178, 178)
bolton.control.cohen

################### ################### ################### ################### ################### 
# Part 3: Some Bayesian work
################### ################### ################### ################### ################### 

################### Update on a prior for mental health interventions derived from the literature ###################

N.samples <- nrow(strongminds.data.1)

# deflate means and widen SDs to attempt to account for attrition
# Note that in practice this actually inflates mean effect sizes, but also blows up SDs
# Uncomment this section and the one following to use this, but I felt inflated means weren't appropriately conservative
#sm.deflations <- deflate.mean.and.sd(strongminds.data.1$pre_mean, strongminds.data.1$pre_sd, strongminds.data.1$post_mean, strongminds.data.1$post_sd, strongminds.data.1$pre_n, strongminds.data.1$post_n)
#strongminds.data.1$post_mean_2 <- sm.deflations$mean
#strongminds.data.1$post_sd_2 <- sm.deflations$sd

# calculate difference between pre and post
#strongminds.data.1$pre.post.adjusted.change <- strongminds.data.1$post_mean_2 - strongminds.data.1$pre_mean
#strongminds.data.1$pre.post.adjusted.sd <- sqrt(strongminds.data.1$pre_sd^2 + strongminds.data.1$post_sd_2^2)

strongminds.data.1$pre.post.adjusted.change <- strongminds.data.1$post_mean - strongminds.data.1$pre_mean
strongminds.data.1$pre.post.adjusted.sd <- sqrt(strongminds.data.1$pre_sd^2 + strongminds.data.1$post_sd^2)

# use the "control effect size" from Bolton 2003, converted into our units
# and adjusted for program length (16 in the case of Bolton, 8 for us)
strongminds.data.1$counterfactual.score.change <- strongminds.data.1$pre.post.adjusted.sd * (bolton.control.cohen/16)*8

# adjusted change (relative to control)
# create backed-out SDs for SM programs taking into account attrition
# assume a null effect for all attriters and recalculate the sd on that basis
strongminds.data.1$adjusted_change <- strongminds.data.1$pre.post.adjusted.change - strongminds.data.1$counterfactual.score.change

# standardized mean difference
strongminds.data.1$standardized.mean.difference <- strongminds.data.1$adjusted_change / strongminds.data.1$pre.post.adjusted.sd

# via https://pubmed.ncbi.nlm.nih.gov/31948935/
# a prior on the standardized mean difference between pre and post in a meta analysis
# of interventions to mitigate depressive symptoms in adults
# 0.94 [0.89, 1.30], so an s.d. of approximately 0.1
# I chose "other psychotherapies" instead of "interpersonal therapy" for adults with depression as the appropriate reference class

prior_mean <- -0.94 # negative because in our context this is a reduction in PHQ-9 scores
prior_sd <- 0.1

sm_data <- list(N = N.samples,
                adjusted_change = strongminds.data.1$standardized.mean.difference,
                #adjusted_change = strongminds.data.1$adjusted_change,
                #change_sd = strongminds.data.1$pre.post.adjusted.sd,
                change_sd = rep(1, length(strongminds.data.1$standardized.mean.difference)),  # because it's standardized
                prior_mean = prior_mean,
                prior_sd = prior_sd
)

fit_rstan <- stan(
  file = "sm_meta1.stan",
  data = sm_data,
  iter = 4000,
  cores = mc.cores
)


extracted.fit <- rstan::extract(fit_rstan)
mu.posterior <- extracted.fit$mu
plot(density(mu.posterior))

# posterior median effect estimate
median(mu.posterior)
# SDs of the score difference are in the neighborhood of 4, so this suggests a reduction of around 4

### Visualize prior, posterior, evidence
prior.density <- rnorm(10000, prior_mean, prior_sd)

ggplot() +
  geom_density(aes(x = prior.density, fill="Prior"), alpha=0.5) +
  geom_density(aes(x = mu.posterior, fill="Posterior"), alpha=0.5) +
  geom_density(aes(x = strongminds.data.1$standardized.mean.difference, fill="Data"), alpha=0.5) +
  guides(fill=guide_legend(title="")) +
  ggtitle("The prior is doing most of the work here!")
  
################### Update on a more conservative prior ###################

conservative_mean <- 0
conservative_sd <- 1

sm_data_con <- list(N = N.samples,
                adjusted_change = strongminds.data.1$standardized.mean.difference,
                #adjusted_change = strongminds.data.1$adjusted_change,
                #change_sd = strongminds.data.1$pre.post.adjusted.sd,
                change_sd = rep(1, length(strongminds.data.1$standardized.mean.difference)),  # because it's standardized
                prior_mean = conservative_mean,
                prior_sd = conservative_sd
)

fit_rstan_con <- stan(
  file = "sm_meta1.stan",
  data = sm_data_con,
  iter = 4000,
  cores = mc.cores
)


extracted.fit.con <- rstan::extract(fit_rstan_con)
mu.posterior.con <- extracted.fit.con$mu
plot(density(mu.posterior.con))

# posterior median effect estimate
median(mu.posterior.con)
# SDs of the score difference are in the neighborhood of 4, so this suggests a reduction of around 7

### Visualize prior, posterior, evidence
prior.density.con <- rnorm(10000, conservative_mean, conservative_sd)

ggplot() +
  geom_density(aes(x = prior.density.con, fill="Prior"), alpha=0.5) +
  geom_density(aes(x = mu.posterior.con, fill="Posterior"), alpha=0.5) +
  geom_density(aes(x = strongminds.data.1$standardized.mean.difference, fill="Data"), alpha=0.5) +
  guides(fill=guide_legend(title="")) +
  ggtitle("The data — which could be hard to believe —\nis doing most of the work here!")



######### UNUSED DETRITUS ########

# use cmdstan if you have to set adapt_delta
# there's bug -- it crashes R with rstan
mod <- cmdstan_model("sm_meta1.stan")
fit_rstan <- mod$sample(
  data = sm_data,
  chains = 4,
  adapt_delta=0.999,
  parallel_chains = mc.cores,
)
