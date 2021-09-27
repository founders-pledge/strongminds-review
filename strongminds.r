# TODO: Cochrane prior
# TODO: Economic growth costs of depression

library(rstan)
library(cmdstanr)
library(tidyverse)

setwd("/users/mattlerner/desktop/evaluations/strongminds")

mc.cores = parallel::detectCores()

################### Read in StrongMinds programmes data ###################

strongminds.data.1 <- read.csv("data/StrongMinds data via HLI - Sheet1.csv")
colnames(strongminds.data.1) <- c('program','type','country','year','pre_n','pre_mean','pre_sd','post_n','post_mean','post_sd','mean_change','change_relative_core')
strongminds.data.1$pre_n <- as.numeric(gsub(",","",strongminds.data.1$pre_n))
strongminds.data.1$post_n <- as.numeric(gsub(",","",strongminds.data.1$post_n))
strongminds.data.1$attrition <- strongminds.data.1$post_n - strongminds.data.1$pre_n
strongminds.data.1$pct_attrition <- strongminds.data.1$attrition/strongminds.data.1$pre_n

# create type indicator
strongminds.data.1 <- strongminds.data.1 %>% mutate(type_id = group_indices(., type))

################### Functions and parameters ###################
N <- 10000

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

# inverse quantile function for exponential distribution
invq.exp <- function(p, quantile) {
  lambda <- (-1 * (log(1-p)))/quantile
  return(lambda)
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

################### Bolton 2003, replicating HLI calculation of Cohen's D ###################

#bolton.cohen.d <- cohen.d(17.47, 3.55, 1.1, 1.1, 107, 117) # whole group
bolton.cohen.d <- cohen.d(11.59, 2.38, 0.8, 0.75, 107, 117) # all eligible persons - this seems to be the HLI number
bolton.cohen.d

# calculate effect size for CONTROL group, e.g. standardized pre-post difference among the untreated
bolton.control.cohen <- cohen.d(23.65, 21.14, 6.3, 8.19, 178, 178)
bolton.control.cohen

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

################### Assessing DALY cost of a reduced PHQ-9 score ###################

average.sm.pretreatment.mean <- mean(strongminds.data.1$pre_mean)

# PDF for distribution of PHQ-9 scores over population, using mean of Zambia & Uganda rate parameters
overall.parameter <- mean(c(zambia.parameter, uganda.parameter))
phq9.dist <- rexp(N, overall.parameter)

# given this distribution, what percentage reduction in major depression symptoms can we expect
# from an expected reduction of 1 point?
mild.depression.cutoff <- 5
moderate.depression.cutoff <- 10
severe.depression.cutoff <- 15 # manual has cutoffs at 15 for "moderately severe" and 20 for "severe". We combine them here.
depression.cutoff <- 10
score.reduction <- 1

# percentage-point incidence in the incidence of major depression
normal.mild.incidence <- mean(phq9.dist >= mild.depression.cutoff & phq9.dist < moderate.depression.cutoff)
normal.moderate.incidence <- mean(phq9.dist >= moderate.depression.cutoff & phq9.dist < severe.depression.cutoff)
normal.severe.incidence <- mean(phq9.dist >= severe.depression.cutoff)

# subtract 1 point from distribution
phq9.modified <- phq9.dist - score.reduction

# modified incidence of each type of depression
modified.mild.incidence <- mean(phq9.modified >= mild.depression.cutoff & phq9.modified < moderate.depression.cutoff)
modified.moderate.incidence <- mean(phq9.modified >= moderate.depression.cutoff & phq9.modified < severe.depression.cutoff)
modified.severe.incidence <- mean(phq9.modified >= severe.depression.cutoff)

mild.change <- normal.mild.incidence - modified.mild.incidence
moderate.change <- normal.moderate.incidence - modified.moderate.incidence
severe.change <- normal.severe.incidence - modified.severe.incidence

mild.change
moderate.change
severe.change


###### StrongMinds phase 2 assessment data
# evaluation
# https://drive.google.com/file/d/1qSQuM_BF8H-ufzfsETyq2c8jYZ-KHo4d/view


###### StrongMinds 2020 RCT data

################### Meta analysis of StrongMinds pre/post data ###################

N <- nrow(strongminds.data.1)

# deflate means and widen SDs to attempt to account for attrition
# TODO: THIS INFLATES EFFECTS!
sm.deflations <- deflate.mean.and.sd(strongminds.data.1$pre_mean, strongminds.data.1$pre_sd, strongminds.data.1$post_mean, strongminds.data.1$post_sd, strongminds.data.1$pre_n, strongminds.data.1$post_n)
strongminds.data.1$post_mean_2 <- sm.deflations$mean
strongminds.data.1$post_sd_2 <- sm.deflations$sd

# calculate difference between pre and post
strongminds.data.1$pre.post.adjusted.change <- strongminds.data.1$post_mean_2 - strongminds.data.1$pre_mean
strongminds.data.1$pre.post.adjusted.sd <- sqrt(strongminds.data.1$pre_sd^2 + strongminds.data.1$post_sd_2^2)


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
# 1.09 [0.89, 1.30], so an s.d. of approximately 0.1
prior_mean <- -1.09
prior_sd <- 0.1
# TODO: the question here is whether we should use a much more skeptical prior

sm_data <- list(N = N,
                adjusted_change = strongminds.data.1$standardized.mean.difference,
                #adjusted_change = strongminds.data.1$adjusted_change,
                #change_sd = strongminds.data.1$pre.post.adjusted.sd,
                change_sd = 1, # because it's standardized
                prior_mean = prior_mean,
                prior_sd = prior_mean
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


######### DETRITUS ########

# use cmdstan if you have to set adapt_delta
# there's bug -- it crashes R with rstan
mod <- cmdstan_model("sm_meta1.stan")
fit_rstan <- mod$sample(
  data = sm_data,
  chains = 4,
  adapt_delta=0.999,
  parallel_chains = mc.cores,
)

