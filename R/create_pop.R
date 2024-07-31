rm(list=ls())
# load packages
library(simPop)
library(tidyverse)

#load and rename data an change age to numeric and pid to factor
data(eusilc13puf, package = "simPop")
??eusilc13puf
data <- eusilc13puf
data <- data |>  mutate(
  age = as.numeric(as.character(age)),
  pid = as.factor(pid)
  ) 

# define input
## approx. 20 seconds computation time
inp <-
  specifyInput(
    data = data,
    hhid = "db030",
    hhsize = "hsize",
    strata = "db040",
    weight = "rb050"
  )
## in the following, nr_cpus are selected automatically

simPop <-
  simStructure(
    data = inp,
    method = "direct",
    basicHHvars = c("age", "rb090", "db040") # Gender & Bundesland
  )
simPop <-
  simCategorical(
    simPop,
    additional = c("pl031", "pb220a", "pb190", "pe040", "pl111"), # economic status & citizenship status
    method = "multinom"#, # also try different methods such as xgboost
   # nr_cpus = 1
  )

# Add known margins
# MT: for simPop we would need to add population margins, not sample margins
#     and: no need for changing the margins for the simulated population.
#     = this is our synthetic truth, and margins are the truth.
# margins <-
#  as.data.frame(xtabs(rep(1, nrow(data)) ~ data$db040 + data$rb090 + data$pb220a))
#
#colnames(margins) <-
#  c("db040", "rb090", "pb220a", "freq") #strata, gender, citizenship status
##margins$freq <- margins$freq
#simPop <- addKnownMargins(simPop, margins)


# define reg model
regModel = ~rb090+hsize+pl031+pb220a+pb190+pe040+pl111 # gender, economic status, citizenship status
# besser?
regModel = ~rb090+hsize+pl031+pb220a+age+pb190+pe040+pl111
# noch besser?
regModel = ~rb090*pb220a+hsize+pl031+age+pb190+pe040+pl111
# multinominal model
simPop <- simContinuous(simPop, additional="pgrossIncome", # used to be netIncome
                         regModel = regModel,
                         upper=200000, equidist=FALSE#, nr_cpus=1
                        ) # also try differnt methods such as xgboost

#regModel = ~rb090+hsize+pl031+pb220a
# better? 
regModel = ~rb090+hsize+pl031+pb220a+pgrossIncome+I(pgrossIncome^2)+pb190+pe040+pl111
# multinominal model
# simPop <- simContinuous(simPop, additional="hgrossIncome", # used to be netIncome
#                        regModel = regModel,
#                        upper=200000, equidist=FALSE#, nr_cpus=1
#                        )

seed<-set.seed(1)

simPop <- simComponents(
  simPop,
  total = "pgrossIncome",
  components = c(
    "py010g",
    "py021g",
    "py050g",
    "py080g",
    "py090g",
    "py100g",
    "py110g",
    "py120g",
    "py130g",
    "py140g"
  ),
  conditional = c("pl031","pb220a","rb090"), # check what is best
  # check what is with getCatName
  replaceEmpty = c("sequential", "min"),
  seed
)


# simPop <- simComponents(
#   simPop,
#   total = "hgrossIncome",
#   components = c(
#     "hy040g",
#     "hy050g",
#     "hy060g",
#     "hy070g",
#     "hy080g",
#     "hy090g",
#     "hy100g",
#     "hy110g",
#     "hy120g",
#     "hy130g",
#     "hy140g"
#   ),
#   conditional = c("pl031", "pb220a","rb090"),
#   # check what is with getCatName
#   replaceEmpty = c("sequential", "min"),
#   seed
# )

# save simPop object
# save(simPop, file = "data/pop.rda")

# save only df
pop <- pop(simPop)
save(pop, file = "data/pop.rda")
# saveRDS(pop, file = "data/pop/Pop_U_star.rds")

# This is a simple approach where the correlation between the phincome and the
# hhincome is lost. For a first analysis this is fine.
# What we can do further is to calculate a mean per household for hgrossincome
# up to row 50 is fine for phh income, do the simulation there, then aggregate 
# eusilc vars on hh level, then simulate this df with (rows 1-50) 
# (weights check, right factor to get the right amount of hh (approx 2-3 mio))
# Then leftjoin hh things, with id



# load packages
library(simPop)
library(simFrame)
library(dplyr)

#load data 
#load("data/pop/Pop_U_star.RData")
load("data/pop.rda")
data(eusilc13puf, package = "simPop")

# setting a seed
set.seed(123)

# Number of households in the original sample:
nhh <- length(unique(eusilc13puf$db030))

# 1. Sample Methods ----

## 1.1 SRSWOR (of households) ----
# s_srswor <- simFrame::draw(simPop@pop@data, 
#                  grouping = "db030",  
#                  size = nhh)
# summary(s_srswor$.weight)
# dim(s_srswor) # ok
# sum(s_srswor$.weight) # ok

## 1.2 Proportional stratified SRSWOR ----
hh <- pop |> 
  group_by(db030) |> 
  summarize(region = unique(db040)) |> 
  group_by(region) |>  
  summarize(counts = n()) |>  
  select(counts)
sizes <- unlist(hh / sum(hh) * nhh)
sum(sizes)

samp <- simFrame::draw(pop,
                                 design = "db040",
                                 grouping = "db030",
                                 size = sizes)
save(samp, file = "data/samp.rda")



