rm(list = ls())
library(simPop)
library(fitdistrplus)
?correctSingleHeap
?correctHeaps
source("~/workspace/heaping/R/correctHeap.R", echo=TRUE)


introduce_heaping <- function(df, age_var, heap_ratio, interval = 5) {
  # Convert to data.table if not already
  is_data_table <- is.data.table(df)
  # if (!is_data_table) {
  #   setDT(df)
  # }
  df <- data.frame(df)
  
  # Check if the age variable exists in the data frame
  if (!age_var %in% names(df)) {
    stop("The specified age variable does not exist in the data frame.")
  }
  
  # Ensure the age variable is numeric
  if (!is.numeric(df[[age_var]])) {
    stop("The specified age variable must be numeric.")
  }
  
  # Determine the number of values to change
  n <- nrow(df)
  n_to_heap <- round(heap_ratio * n)
  
  # Sample indices to heap
  indices_to_heap <- sample(n, n_to_heap)
  
  # Round ages to the nearest specified interval
  #df[indices_to_heap, (age_var) := round(df[[age_var]][indices_to_heap] / interval) * interval]
  df[indices_to_heap, age_var] <- round(df[indices_to_heap, age_var] / interval) * interval
  
  # Convert back to data.frame if necessary
  # if (!is_data_table) {
  #   setDF(df)
  # }
  if (!is_data_table) {
    df <- data.table(df)
  }
  return(df)
}



## create some artificial data
age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
age <- round(age[age < 93])
barplot(table(age))

## artificially introduce age heaping and correct it:
# heaps every 5 years
year5 <- seq(0, max(age), 5)
age5 <- sample(c(age, age[age %in% year5]))
cc5 <- rep("darkgrey", length(unique(age)))
cc5[year5+1] <- "yellow"

pdf("figures/age5heaps1.pdf")
barplot(table(age5), col=cc5)
dev.off()
pdf("figures/age5heaps2.pdf")
barplot(table(correctHeaps2(age5, heaps="5year", method="lnorm")), col=cc5)
dev.off()


# heaps every 10 years
year10 <- seq(0, max(age), 10)
age10 <- sample(c(age, age[age %in% year10]))
cc10 <- rep("darkgrey", length(unique(age)))
cc10[year10+1] <- "yellow"
pdf("figures/age10heaps1.pdf")
barplot(table(age10), col=cc10)
dev.off()

pdf("figures/age10heaps2.pdf")
barplot(table(correctHeaps(age10, heaps="10year", method="lnorm")), col=cc10)
dev.off()

# the first 5 observations should be unchanged
pdf("figures/age10heaps2first5.pdf")
barplot(table(correctHeaps(age10, heaps="10year", method="lnorm", fixed=1:5)), col=cc10)
dev.off()

## Single heap

## create some artificial data
age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
age <- round(age[age < 93])
barplot(table(age))

## artificially introduce an age heap for a specific year
## and correct it
age23 <- c(age, rep(23, length=sum(age==23)))
cc23 <- rep("darkgrey", length(unique(age)))
cc23[24] <- "yellow"
pdf("figures/age1heap1.pdf")
barplot(table(age23), col=cc23)
dev.off()
pdf("figures/age1heap1.pdf")
barplot(table(correctSingleHeap(age23, heap=23, before=2, after=3, method="lnorm")), col=cc23)
dev.off()
pdf("figures/age1heap2.pdf")
barplot(table(correctSingleHeap(age23, heap=23, before=5, after=5, method="lnorm")), col=cc23)
dev.off()
pdf("figures/age1heap2first5.pdf")
# the first 5 observations should be unchanged
barplot(table(correctSingleHeap(age23, heap=23, before=5, after=5, method="lnorm",
                                fixed=1:5)), col=cc23)
dev.off()

##########################################
###.  Chicago.  ##########################
##########################################

chic <- function(chicago, form = formula("ln.real.wage ~ age + female + LTHS +
                   some.college + college + advanced.degree | foreign.born | LTHS +
                   some.college + college + advanced.degree"), kind = "barplot",
                 seed = 911){
  
  if(!is.null(seed)) set.seed(seed)
  library(data.table)
  chicago_heap <- chicago
  chicago_heap <- introduce_heaping(chicago, "age", 0.27, interval = 5)
  chicago_after_model <- chicago_after <- chicago_heap
  chicago_after$age <- correctHeaps2(chicago_heap$age, heaps = "5year", method = "lnorm")
  chicago_after_model$age <- correctHeaps2(chicago_heap$age, heaps = "5year", method = "lnorm",
                                           model = as.formula("age ~ ln.real.wage + female + 
                                                    some.college + college + advanced.degree + foreign.born + LTHS +
                    some.college + college + advanced.degree"), dataModel = chicago_heap)
  
  res <- NULL
  
  if(kind == "barplot"){
    t1 <- table(chicago$age)
    t2 <- table(chicago_heap$age)
    t3 <- table(chicago_after$age)
    t4 <- table(chicago_after_model$age)
    xlims <- max(c(t1, t2, t3, t4)) 
    par(mfrow = c(2,2))
    barplot(t1, ylim = c(0, xlims), xlab = "age (original)")
    barplot(t2, ylim = c(0, xlims), xlab = "age (heaping)")
    barplot(t3, ylim = c(0, xlims), xlab = "age (corrected)")
    barplot(t4, ylim = c(0, xlims), xlab = "age (corrected + model)")
  }
  
  if(kind == "densityplot"){
    d <- density(chicago$age, bw = 1)
    par(mar = c(4,4,0.3,0.1))
    plot(d, main = "", xlab = "age", lwd = 2)
    lines(density(chicago_heap$age, bw = d$bw), col = "red")
    lines(density(chicago_after$age, bw = d$bw), col = "blue")
    lines(density(chicago_after_model$age, bw = d$bw), col = "magenta")
    legend("topright", legend = c("original", "heaping", "corrected", "corrected (model)"),
           col = c("black", "red", "blue", "magenta"), lwd = c(2,1,1,1))
  }
  
  if(kind %in% c("twofold", "beta.diff")){
    results <- oaxaca(formula = form, data = chicago, R = 1000)
    results_heap <- oaxaca(formula = form, data = chicago_heap, R = 1000)
    results_after <- oaxaca(formula = form, data = chicago_after, R = 1000)
    results_after_model <- oaxaca(formula = form, data = chicago_after_model, R = 1000)
  }
  
  if(kind == "twofold"){
    variables <- c("age", "female", "college")
    columns <- c("group.weight", "coef(unexplained A)", "coef(unexplained B)")
    res <- list("chicago" = results$twofold$variables[[5]][variables, columns],
      "chicago_heap" = results_heap$twofold$variables[[5]][variables, columns],
      "chicago_after" = results_after$twofold$variables[[5]][variables, columns],
      "chicago_after_model" = results_after_model$twofold$variables[[5]][variables, columns])
  }
  
  if(kind == "beta.diff"){
   res <- c(results$beta$beta.diff["age"],
    results_heap$beta$beta.diff["age"],
    results_after$beta$beta.diff["age"],
    results_after_model$beta$beta.diff["age"])
   names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  
  if(kind %in% c("coef", "conf.int")){
    lm1 <- lm(formula = ln.real.wage ~ age + female + LTHS +
                some.college + college + advanced.degree * foreign.born * LTHS +
                some.college + college + advanced.degree, data = chicago)
#    coef(lm1)[2]
    
    lm1_heap <- lm(formula = ln.real.wage ~ age + female + LTHS +
                     some.college + college + advanced.degree * foreign.born * LTHS +
                     some.college + college + advanced.degree, data = chicago_heap)

    lm1_after <- lm(formula = ln.real.wage ~ age + female + LTHS +
                      some.college + college + advanced.degree * foreign.born * LTHS +
                      some.college + college + advanced.degree, data = chicago_after)

    lm1_after_model <- lm(formula = ln.real.wage ~ age + female + LTHS +
                            some.college + college + advanced.degree * foreign.born * LTHS +
                            some.college + college + advanced.degree, data = chicago_after_model)
  }
  
  if(kind == "coef"){
    res <- c(coef(lm1)[2], coef(lm1_heap)[2], coef(lm1_after)[2], coef(lm1_after_model)[2])
    names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  
  if(kind == "conf.int"){
    res <- list("chicago" = confint(lm1)[2, ],
                          "chicago_heap" = confint(lm1_heap)[2, ],
                          "chicago_after" = confint(lm1_after)[2, ],
                          "chicago_after_model" = confint(lm1_after_model)[2, ])
    res$kind <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }

  if(kind %in% c("RMSE", "Rsquared", "MAE", "RMSESD", "MAESD")){
    require(caret)
    chicago$agecut <- as.numeric(as.character(cut(chicago$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    chicago_heap$agecut <- as.numeric(as.character(cut(chicago_heap$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    chicago_after$agecut <- as.numeric(as.character(cut(chicago_after$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    chicago_after_model$agecut <- as.numeric(as.character(cut(chicago_after_model$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    # Define the formula for the model
    formula <- agecut ~ ln.real.wage + 
      female + 
      LTHS * foreign.born * advanced.degree + 
      some.college + 
      college
    # Define the training control with 10-fold cross-validation and 3 repetitions
    train_control <- trainControl(method = "repeatedcv", 
                                  number = 10, 
                                  repeats = 3)
    lm_chicago_age_caret <- train(formula, 
                                  data = na.omit(chicago), 
                                  method = "lm", 
                                  trControl = train_control)
    lm_chicago_heap_age_caret <- train(formula, 
                                       data = na.omit(chicago_heap), 
                                       method = "lm", 
                                       trControl = train_control)
    lm_chicago_after_age_caret <- train(formula, 
                                        data = na.omit(chicago_after), 
                                        method = "lm", 
                                        trControl = train_control)
    lm_chicago_after_model_age_caret <- train(formula, 
                                              data = na.omit(chicago_after_model), 
                                              method = "lm", 
                                              trControl = train_control)

  }
  
  if(kind == "RMSE"){
    res <- c(lm_chicago_age_caret$results$RMSE, 
             lm_chicago_heap_age_caret$results$RMSE, 
             lm_chicago_after_age_caret$results$RMSE, 
             lm_chicago_after_model_age_caret$results$RMSE)
    names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  if(kind == "Rsquared"){
    res <- c(lm_chicago_age_caret$results$Rsquared, 
             lm_chicago_heap_age_caret$results$Rsquared, 
             lm_chicago_after_age_caret$results$Rsquared, 
             lm_chicago_after_model_age_caret$results$Rsquared)
    names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  if(kind == "MAE"){
    res <- c(lm_chicago_age_caret$results$MAE, 
             lm_chicago_heap_age_caret$results$MAE, 
             lm_chicago_after_age_caret$results$MAE, 
             lm_chicago_after_model_age_caret$results$MAE)
    names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  if(kind == "RMSESD"){
    res <- c(lm_chicago_age_caret$results$RMSESD, 
             lm_chicago_heap_age_caret$results$RMSESD, 
             lm_chicago_after_age_caret$results$RMSESD, 
             lm_chicago_after_model_age_caret$results$RMSESD)
    names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  if(kind == "MAESD"){
    res <- c(lm_chicago_age_caret$results$MAESD, 
             lm_chicago_heap_age_caret$results$MAESD, 
             lm_chicago_after_age_caret$results$MAESD, 
             lm_chicago_after_model_age_caret$results$MAESD)
    names(res) <- c("chicago", "chicago_heap", "chicago_after", "chicago_after_model")
  }
  return(res)
}


chic(chicago)
dev.off()
chic(chicago, kind = "densityplot")
chic(chicago, kind = "twofold")
chic(chicago, kind = "beta.diff")
chic(chicago, kind = "coef")
chic(chicago, kind = "conf.int")
chic(chicago, kind = "RMSE")
chic(chicago, kind = "Rsquared")
chic(chicago, kind = "MAE")
chic(chicago, kind = "RMSESD")
chic(chicago, kind = "MAESD")






###  Pop --> Sample --> correct variances.

library(laeken)
data(ses)

# Oaxana-Blinder decomposition

library(oaxaca)
data("chicago")
results <- oaxaca(formula = ln.real.wage ~ age + female + LTHS +
                    some.college + college + advanced.degree | foreign.born | LTHS +
                    some.college + college + advanced.degree, data = chicago, R = 1000)
results
#plot(results, components = c("endowments","coefficients"))
summary(results$reg$reg.pooled.2)$coefficients["LTHS",]
results$beta$beta.diff["age"]
#plot(results, decomposition = "twofold", group.weight = -1)
#plot(results, decomposition = "twofold", weight = -1,
#        unexplained.split = TRUE, components = c("unexplained A",
#                                                   "unexplained B"), component.labels = c("unexplained A" =
#                                                                                              "In Favor of Natives", "unexplained B" = "Against the Foreign-Born"),
#        variables = c("age", "female", "college"), variable.labels = c("age" =
#                                                                           "Years of Age", "female" = "Female", "college" = "College Education"),
#     group.weight = -1)

#plot(results, decomposition = "twofold", weight = -1,
#      unexplained.split = TRUE, components = c("unexplained A",
#                                                 "unexplained B"), component.labels = c("unexplained A" =
#                                                                                            "In Favor of Natives", "unexplained B" = "Against the Foreign-Born"),
#      component.left = TRUE, variables = c("age","female","college"),
#      variable.labels = c("age" = "Years of Age", "female" = "Female",
#                            "college" = "College Education"),group.weight = -1)

variables <- c("age", "female", "college")
columns <- c("group.weight", "coef(unexplained A)", "coef(unexplained B)")
results$twofold$variables[[5]][variables, columns]

## now with heaping
library(data.table)
chicago_heap <- chicago
chicago_heap <- introduce_heaping(chicago, "age", 0.27, interval = 5)
results_heap <- oaxaca(formula = ln.real.wage ~ age + female + LTHS +
                    some.college + college + advanced.degree | foreign.born | LTHS +
                    some.college + college + advanced.degree, data = chicago_heap, R = 1000)
#plot(results_heap, components = c("endowments","coefficients"))
results$twofold$variables[[5]][variables, columns]
results_heap$twofold$variables[[5]][variables, columns]

chicago_after_model <- chicago_after <- chicago_heap
source("~/workspace/heaping/R/correctHeap.R", echo=TRUE)
chicago_after$age <- correctHeaps2(chicago_heap$age, heaps = "5year", method = "lnorm")


## test model

chicago_after_model$age <- correctHeaps2(chicago_heap$age, heaps = "5year", method = "lnorm",
                                   model = as.formula("age ~ ln.real.wage + female + 
                                                    some.college + college + advanced.degree + foreign.born + LTHS +
                    some.college + college + advanced.degree"), dataModel = chicago_heap)

##

results_after <- oaxaca(formula = ln.real.wage ~ age + female + LTHS +
                         some.college + college + advanced.degree | foreign.born | LTHS +
                         some.college + college + advanced.degree, data = chicago_after, R = 1000)
results_after_model <- oaxaca(formula = ln.real.wage ~ age + female + LTHS +
                          some.college + college + advanced.degree | foreign.born | LTHS +
                          some.college + college + advanced.degree, data = chicago_after_model, R = 1000)
#plot(results_after, components = c("endowments","coefficients"))
results$twofold$variables[[5]][variables, columns]
results_heap$twofold$variables[[5]][variables, columns]
results_after$twofold$variables[[5]][variables, columns]
results_after_model$twofold$variables[[5]][variables, columns]
results$beta$beta.diff["age"]
results_heap$beta$beta.diff["age"]
results_after$beta$beta.diff["age"]
results_after_model$beta$beta.diff["age"]
plot(density(chicago$age))
lines(density(chicago_heap$age), col = "red")
lines(density(chicago_after$age), col = "blue")
lines(density(chicago_after_model$age), col = "magenta")


par(mfrow = c(2,2))
barplot(table(chicago$age))
barplot(table(chicago_heap$age))
barplot(table(chicago_after$age))
barplot(table(chicago_after_model$age))

lm1 <- lm(formula = ln.real.wage ~ age + female + LTHS +
                         some.college + college + advanced.degree * foreign.born * LTHS +
                         some.college + college + advanced.degree, data = chicago)
summary(lm1)
coef(lm1)[2]
confint(lm1)[2, ]
lm1_heap <- lm(formula = ln.real.wage ~ age + female + LTHS +
                 some.college + college + advanced.degree * foreign.born * LTHS +
                 some.college + college + advanced.degree, data = chicago_heap)
coef(lm1_heap)[2]
confint(lm1_heap)[2, ]
lm1_after <- lm(formula = ln.real.wage ~ age + female + LTHS +
                 some.college + college + advanced.degree * foreign.born * LTHS +
                 some.college + college + advanced.degree, data = chicago_after)
coef(lm1_after)[2]
confint(lm1_after)[2, ]
lm1_after_model <- lm(formula = ln.real.wage ~ age + female + LTHS +
                  some.college + college + advanced.degree * foreign.born * LTHS +
                  some.college + college + advanced.degree, data = chicago_after_model)
coef(lm1_after_model)[2]
confint(lm1_after_model)[2, ]

## regression error on age:


library(caret)

chicago$agecut <- as.numeric(as.character(cut(chicago$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
chicago_heap$agecut <- as.numeric(as.character(cut(chicago_heap$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
chicago_after$agecut <- as.numeric(as.character(cut(chicago_after$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
chicago_after_model$agecut <- as.numeric(as.character(cut(chicago_after_model$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))


# Define the formula for the model
formula <- agecut ~ ln.real.wage + 
  female + 
  LTHS * foreign.born * advanced.degree + 
  some.college + 
  college

# Define the training control with 10-fold cross-validation and 3 repetitions
train_control <- trainControl(method = "repeatedcv", 
                              number = 10, 
                              repeats = 3)

# Fit the model using the train function from caret
# lm_pop_age <- train(formula, 
#                     data = samp, 
#                     method = "lm", 
#                     trControl = train_control)
lm_chicago_age_caret <- train(formula, 
                           data = na.omit(chicago), 
                           method = "lm", 
                           trControl = train_control)
lm_chicago_heap_age_caret <- train(formula, 
                                data = na.omit(chicago_heap), 
                                method = "lm", 
                                trControl = train_control)
lm_chicago_after_age_caret <- train(formula, 
                                   data = na.omit(chicago_after), 
                                   method = "lm", 
                                   trControl = train_control)
lm_chicago_after_model_age_caret <- train(formula, 
                                    data = na.omit(chicago_after_model), 
                                    method = "lm", 
                                    trControl = train_control)
# Print the results
print(lm_chicago_age_caret)
print(lm_chicago_heap_age_caret)
print(lm_chicago_after_age_caret)
print(lm_chicago_after_model_age_caret)

## train on sample, evaluate on pop

pred_noheap <- predict(lm_samp_age, newdata = pop)
pred_heap <- predict(lm_samp_heap_age, newdata = pop)

sqrt(mean((pop$age - pred_noheap)^2))
sqrt(mean((pop$age - pred_heap)^2))
head(pred_noheap)
head(pred_heap)
summary(lm_samp_age)
summary(lm_samp_heap_age)
plot(density(pred_noheap))
lines(density(pred_heap))


regress <- lm(log(earningsHour) ~ age + 
              education + occupation + contract + fullPart + lengthService + NACE1 + location + payAgreement + size + economicFinanc, 
            weight = weights, 
            data = ses)
summary(regress)

## SILC

library(dplyr)
load("data/pop.rda")
load("data/samp.rda")
# data("eusilc13puf", package = "simPop")

# Define the function
handle_missing_values <- function(df) {
  df <- df %>%
    mutate(across(where(is.factor), ~ ifelse(is.na(.), "not applicable", as.character(.))),
           across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))
  
  # Convert modified factor columns back to factor type
  df <- df %>%
    mutate(across(where(is.character), ~ as.factor(.)))
  
  return(df)
}

pop <- handle_missing_values(pop)
samp <- handle_missing_values(samp)
samp$age <- as.numeric(as.character(samp$age))
pop$age <- as.numeric(as.character(pop$age))

pop <- pop[pop$age >= 18 & pop$age <= 65, ]
samp <- samp[samp$age >= 18 & samp$age <= 65, ]
lm_pop <- lm(log(pgrossIncome+1) ~ age +
                rb090 + 
                pl031 + # Derzeitige Hauptaktivität
                pb220a + 
                pb190 + # Familienstand + 
                pe040 + # Höchster Bildungsabschluss 
                pl111, #+ # Letzter Wirtschaftszweig/Branche (ehemals Erwerbstätige) NACE 
              #     p118000i, # Höchster Bildungsabschluss
              data = pop)
summary(lm_pop)
confint(lm_pop)[2, ]
lm_samp <- lm(log(pgrossIncome+1) ~ age +
             rb090 + 
             pl031 + # Derzeitige Hauptaktivität
             pb220a + 
             pb190 + # Familienstand + 
             pe040 + # Höchster Bildungsabschluss 
             pl111, #+ # Letzter Wirtschaftszweig/Branche (ehemals Erwerbstätige) NACE 
        #     p118000i, # Höchster Bildungsabschluss
          data = samp)
summary(lm_samp)
confint(lm_samp)[2, ]

# Now introduce heaping in age
library(data.table)




samp_heap <- samp
set.seed(123)
samp_heap <- introduce_heaping(samp, "age", 0.27, interval = 10)

lm_samp_heap <- lm(log(pgrossIncome+1) ~ age +
                rb090 + 
                pl031 + # Derzeitige Hauptaktivität
                pb220a + 
                pb190 + # Familienstand + 
                pe040 + # Höchster Bildungsabschluss 
                pl111, #+ # Letzter Wirtschaftszweig/Branche (ehemals Erwerbstätige) NACE 
              #     p118000i, # Höchster Bildungsabschluss
              data = samp_heap)
summary(lm_samp_heap)
confint(lm_samp_heap)[2, ]


## Predict age
# lm_pop_age <- lm(age ~ log(pgrossIncome + 1) +
#                      rb090 + 
#                      pl031 + 
#                      pb220a + 
#                      pb190 + 
#                      pe040 + 
#                      pl111,
#                    data = pop)
lm_samp_age <- lm(age ~ log(pgrossIncome + 1) +
                     rb090 + 
                     pl031 + # Derzeitige Hauptaktivität
                     pb220a + 
                     pb190 + # Familienstand + 
                     pe040 + # Höchster Bildungsabschluss 
                     pl111, #+ # Letzter Wirtschaftszweig/Branche (ehemals Erwerbstätige) NACE 
                   #     p118000i, # Höchster Bildungsabschluss
                   data = samp)
lm_samp_heap_age <- lm(age ~ log(pgrossIncome + 1) +
                     rb090 + 
                     pl031 + # Derzeitige Hauptaktivität
                     pb220a + 
                     pb190 + # Familienstand + 
                     pe040 + # Höchster Bildungsabschluss 
                     pl111, #+ # Letzter Wirtschaftszweig/Branche (ehemals Erwerbstätige) NACE 
                   #     p118000i, # Höchster Bildungsabschluss
                   data = samp_heap)
summary(lm_samp_age)
summary(lm_samp_heap_age)

library(caret)

# Define the formula for the model
formula <- age ~ log(pgrossIncome + 1) + 
  rb090 + 
  pl031 + 
  pb220a + 
  pb190 + 
  pe040 + 
  pl111

# Define the training control with 10-fold cross-validation and 3 repetitions
train_control <- trainControl(method = "repeatedcv", 
                              number = 10, 
                              repeats = 3)

# Fit the model using the train function from caret
# lm_pop_age <- train(formula, 
#                     data = samp, 
#                     method = "lm", 
#                     trControl = train_control)
lm_samp_age_caret <- train(formula, 
                    data = samp, 
                    method = "lm", 
                    trControl = train_control)
lm_samp_heap_age_caret <- train(formula, 
                     data = samp_heap, 
                     method = "lm", 
                     trControl = train_control)

# Print the results
print(lm_samp_age_caret)
print(lm_samp_heap_age_caret)


## train on sample, evaluate on pop

pred_noheap <- predict(lm_samp_age, newdata = pop)
pred_heap <- predict(lm_samp_heap_age, newdata = pop)

sqrt(mean((pop$age - pred_noheap)^2))
sqrt(mean((pop$age - pred_heap)^2))
head(pred_noheap)
head(pred_heap)
summary(lm_samp_age)
summary(lm_samp_heap_age)
plot(density(pred_noheap))
lines(density(pred_heap))
