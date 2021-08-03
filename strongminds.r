library(rstan)
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

invq.exp <- function(p, quantile) {
  lambda <- (-1 * (log(1-p)))/quantile
  return(lambda)
}

################### Bolton 2003, replicating HLI calculation of Cohen's D ###################

#bolton.cohen.d <- cohen.d(17.47, 3.55, 1.1, 1.1, 107, 117) # whole group
bolton.cohen.d <- cohen.d(11.59, 2.38, 0.8, 0.75, 107, 117) # all eligible persons - this seems to be the HLI number
bolton.cohen.d

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

# percentage-point difference in the incidence of major depression
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

################### Meta analysis of StrongMinds pre/post data ###################

# TODO HERE: attempt to correct for selection effects and regression to the mean

# this should come from bolton 2003
estimated_control_effect <- 2.52
N <- nrow(strongminds.data.1)

sm_data <- list(N = N,
                N_types = max(strongminds.data.1$type_id),
                types = strongminds.data.1$type_id,
                mean_change = strongminds.data.1$mean_change,
                pct_attrition = strongminds.data.1$pct_attrition,
                pre_n = strongminds.data.1$pre_n,
                estimated_control_effect = rep(estimated_control_effect, N))

fit_rstan <- stan(
  file = "sm_meta1.stan",
  data = sm_data
)

###### StrongMinds phase 2 assessment data

# TKTKTK

###### StrongMinds 2020 RCT data