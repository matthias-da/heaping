rm(list = ls())
library(simPop)
library(fitdistrplus)
library(data.table)
library(reshape2)
library(ggplot2)
?correctSingleHeap
?correctHeaps
source("~/workspace/heaping/R/correctHeap.R", echo=TRUE)


source("~/workspace/heaping/R/cmd-functions.R", echo=TRUE)

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



genDat <- function(ratio = 0.27,strangeness = 0.01){
  ## create some artificial data
  age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
  age <- round(age[age < 93])
  simulated_data <- simulate_covariates(age, strangeness)
  simulated_data$income <- log(simulated_data$income)
  
  ## artificially introduce age heaping and correct it:
  # heaps every 5 years
  year5 <- seq(0, max(age), 5)
  age5 <- sample(c(age, age[age %in% year5]), size = length(age))
#  cc5 <- rep("darkgrey", length(unique(age)))
#  cc5[year5+1] <- "yellow"
#  # heaps every 10 years
#  year10 <- seq(0, max(age), 10)
#  age10 <- sample(c(age, age[age %in% year10]))
#  cc10 <- rep("darkgrey", length(unique(age)))
#  cc10[year10+1] <- "yellow"
  age5_heaped <- introduce_heaping(data.frame("age" = age), "age", heap_ratio = ratio, interval = 5)
  age5_corrected <- correctHeaps(as.numeric(unlist(age5_heaped)))
  age5_corrected2 <- correctHeaps2(as.numeric(unlist(age5_heaped)), model = formula("age ~ income + marital_status + education"),
                                   dataModel = simulated_data)

  library("mgcv")
  df <- data.frame(x = 0:95, y = as.numeric(fill_missing_ages(table(age5_heaped))))
  m_gam <- gam(y ~ s(x), data = df)
  p_gam <- round(predict(m_gam, newdata = df, type = "response"))
  
  
  
  return(list("age" = age,
              "age_heaped" = as.numeric(unlist(age5_heaped)),
              "age_corrected" = age5_corrected,
              "age_corrected_model" = age5_corrected2, 
              "age_gam" = p_gam))
}

simSimple <- function(ratio = 0.27, strangeness = 0.01){
  x <- genDat(ratio, strangeness)
  # JSD
  # Convert counts to probability distributions
  tage <- fill_missing_ages(table(x$age))
  tage5 <- fill_missing_ages(table(x$age_heaped))
  tage5_corrected <- fill_missing_ages(table(x$age_corrected))
  tage5_corrected_model <- fill_missing_ages(table(x$age_corrected_model))
  tage5_gam <- x$age_gam
  p1 <- (tage) / sum((tage))
  p2 <- (tage5) / sum((tage5))
  p3 <- (tage5_corrected) / sum((tage5_corrected))
  p3m <- (tage5_corrected_model) / sum((tage5_corrected_model))
  p4 <- tage5_gam / sum(tage5_gam)
  # Compute Jensen-Shannon Divergence
  KLDiv <- function(A, B, epsilon = 1e-10) {
    A <- A + epsilon
    B <- B + epsilon
    sum(A * log(A / B))
  }
  KLDiv_bayes <- function(A, B, epsilon = 1e-10) {
    A <- A + epsilon
    B <- B + epsilon
    length(A) / 2 * log(mean(A / B, na.rm = TRUE) * mean(B / A, na.rm = TRUE))
  }
  JSDiv <- function(A, B) {
    M <- (A + B) / 2
    jsd <- 0.5 * KLDiv(A, M) + 0.5 * KLDiv(B, M)
    return(jsd)
  }
  JSDiv_bayes <- 
  function(A, B) {
    M <- (A + B) / 2
    jsd <- 0.5 * KLDiv_bayes(A, M) + 0.5 * KLDiv_bayes(B, M)
    return(jsd)
  }
  jsd5 <- JSDiv(p1, p2)
  jsd5corrected <- JSDiv(p1, p3)
  jsd5corrected_model <- JSDiv(p1, p3m)
  jsd5gam <- JSDiv(p1, p4)
  jsd_bayes <- JSDiv_bayes(p1, p2)
  jsd5corrected_bayes <- JSDiv_bayes(p1, p3) 
  jsd5corrected_model_bayes <- JSDiv_bayes(p1, p3m) 
  jsd5gam_bayes <- JSDiv_bayes(p1, p4)   
  kld <- KLDiv(p1, p2)
  kld5corrected <- KLDiv(p1, p3)
  kld5corrected_model <- KLDiv(p1, p3m)
  kld5gam <- KLDiv(p1, p4)
  kld_bayes <- KLDiv_bayes(p1, p2)
  kld5corrected_bayes <- KLDiv_bayes(p1, p3)
  kld5corrected_model_bayes <- KLDiv_bayes(p1, p3m)
  kld5gam_bayes <- KLDiv_bayes(p1, p4)
  
  # Compute Mean Absolute Error
  mae <- mean(abs(tage - tage5))
  mae5corrected <- mean(abs(tage - tage5_corrected))
  mae5corrected_model <- mean(abs(tage - tage5_corrected_model))
  mae5gam <- mean(abs(tage - tage5_gam))
  # Calculate the absolute percentage errors
  tage <- tage + 1
  tage5 <- tage5 + 1
  tage5_corrected <- tage5_corrected + 1
  tage5_corrected_model <- tage5_corrected_model + 1
  tage5_gam <- tage5_gam + 1
  ape <- (abs((tage) - (tage5)) / (tage)) * 100
  ape5corrected <- (abs((tage) - (tage5_corrected)) / (tage)) * 100
  ape5corrected_model <- (abs((tage) - (tage5_corrected_model)) / (tage)) * 100
  ape5gam <- (abs((tage) - (tage5_gam)) / (tage)) * 100
  mape <- mean(ape)
  mape5_corrected <- mean(ape5corrected)
  mape5_corrected_model <- mean(ape5corrected_model)
  mape5_gam <- mean(ape5gam)  
  # Perform Chi-Square test
  suppressWarnings(chisq_age5 <- chisq.test(tage, tage5))
  suppressWarnings(chisq_age5_corrected <- chisq.test(tage, tage5_corrected))
  suppressWarnings(chisq_age5_corrected_model <- chisq.test(tage, tage5_corrected_model))
  suppressWarnings(chisq_age5_gam <- chisq.test(tage, tage5_gam))
  
  w <- whipple(x$age)
  w_heaped <- whipple(as.numeric(unlist(x$age_heaped)))
  w_corrected <- whipple(x$age_corrected)
  w_corrected_model <- whipple(x$age_corrected_model)
  w_gam <- whipple(x$age_gam)
  
  res <- c(mae, mae5corrected, mae5corrected_model, mae5gam,
           mape, mape5_corrected, mape5_corrected_model, mape5_gam,
           chisq_age5$statistic, chisq_age5_corrected$statistic, chisq_age5_corrected_model$statistic, chisq_age5_gam$statistic,
           jsd5, jsd5corrected, jsd5corrected_model, jsd5gam,
           jsd_bayes, jsd5corrected_bayes, jsd5corrected_model_bayes, jsd5gam_bayes, 
           kld, kld5corrected, kld5corrected_model, kld5gam,
           kld_bayes, kld5corrected_bayes, kld5corrected_model_bayes, kld5gam_bayes,
           w, w_heaped, w_corrected, w_corrected_model, w_gam)
  names(res) <- c("mae_heaped", "mae_corrected", "mae_corrected_model", "mae_gam",
                  "mape_heaped", "mape_corrected", "mape_corrected_model", "mape_gam",
                  "chisq heaped", "chisq corrected", "chisq corrected_model", "chisq gam",
                  "Jensen-Shannon divergence heaped", "Jensen-Shannon divergence corrected","Jensen-Shannon divergence corrected_model","Jensen-Shannon divergence gam",
                  "Jensen-Shannon divergence Bayes heaped", "Jensen-Shannon divergence Bayes corrected","Jensen-Shannon divergence Bayes corrected_model","Jensen-Shannon divergence Bayes gam",
                  "Kullback-Leibler heaped", "Kullback-Leibler corrected","Kullback-Leibler corrected_model","Kullback-Leibler gam",
                  "Kullback-Leibler Bayes heaped", "Kullback-Leibler Bayes corrected", "Kullback-Leibler Bayes corrected_model","Kullback-Leibler Bayes gam",
                  "Whipple", "Whipple heaped", "Whipple corrected", "Whipple corrected_model", "Whipple gam")
  res
}

simSimple()
set.seed(11)
r <- replicate(10, simSimple())

par(mar = c(4,7.25,0.1,0.1))
mapes <- t(r)[, c("mape_heaped", "mape_corrected", "mape_corrected_model", "mape_gam")]
colnames(mapes) <- c("heaped", "corrected", "corrected (model)", "smoothed (gam)")
pdf("figures/mapes027.pdf")
boxplot(mapes, horizontal = TRUE, las = 1, xlab = "MAPE")
dev.off()

rl <- list()
j <- 0
s <- seq(0, 0.35, 0.025)
for(i in s){
  j <- j + 1
  rl[[j]] <- as.data.frame(t(rowMeans(replicate(100, simSimple(ratio = i, strangeness = 0.1)))))
}
names(rl) <- as.character(s)
df <- rbindlist(rl, use.names = TRUE, fill = TRUE, idcol = "heaping ratio")
df_mape <- df[, c("mape_heaped", "mape_corrected", "mape_corrected_model", "mape_gam", "heaping ratio")]
colnames(df_mape) <- c("heaped", "corrected", "corrected (model)", "smoothed (gam)", "heaping ratio")
df_mape <- reshape2::melt(df_mape, id.vars = "heaping ratio")
df_mape
df_mape$`heaping ratio` <- as.numeric(df_mape$`heaping ratio`)
df_mape$variable <- as.factor(df_mape$variable)

pdf("figures/mapesim.pdf")
ggplot(df_mape, aes(x = `heaping ratio`, y = value, color = variable)) +
  geom_line() +
  geom_point() +  # Add points for better visibility
  labs(#title = "MAPE by Heaping Ratio and Correction Method",
       x = "Heaping Ratio",
       y = "MAPE Value") +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()


df_jsd <- df[, c("Jensen-Shannon divergence heaped", "Jensen-Shannon divergence corrected", "Jensen-Shannon divergence corrected_model", "Jensen-Shannon divergence gam", "heaping ratio")]
colnames(df_jsd) <- c("heaped", "corrected", "corrected (model)","smoothed (gam)", "heaping ratio")
df_jsd <- reshape2::melt(df_jsd, id.vars = "heaping ratio")
df_jsd
df_jsd$`heaping ratio` <- as.numeric(df_jsd$`heaping ratio`)
df_jsd$variable <- as.factor(df_jsd$variable)

pdf("figures/jsdsim.pdf")
ggplot(df_jsd, aes(x = `heaping ratio`, y = value, color = variable)) +
  geom_line() +
  geom_point() +  # Add points for better visibility
  labs(#title = "MAPE by Heaping Ratio and Correction Method",
    x = "Heaping Ratio",
    y = "Jensen-Shannon divergence") +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()

df_whipple <- df[, c("Whipple", "Whipple heaped", "Whipple corrected", "Whipple corrected_model","Whipple gam", "heaping ratio")]
colnames(df_whipple) <- c("non-heaped","heaped", "corrected",  "corrected (model)", "smoothed (gam)", "heaping ratio")
df_whipple <- reshape2::melt(df_whipple, id.vars = "heaping ratio")
df_whipple
df_whipple$`heaping ratio` <- as.numeric(df_whipple$`heaping ratio`)
df_whipple$variable <- as.factor(df_whipple$variable)

pdf("figures/whipplesim.pdf")
ggplot(df_whipple, aes(x = `heaping ratio`, y = value, color = variable)) +
  geom_line() +
  geom_point() +  # Add points for better visibility
  labs(#title = "MAPE by Heaping Ratio and Correction Method",
    x = "Heaping Ratio",
    y = "Whipple index") +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()




##########################################
###.  Chicago.  ##########################
##########################################

library(oaxaca)
data("chicago")
chic <- function(df, 
                 form = formula("ln.real.wage ~ age + female + LTHS +
                   some.college + college + advanced.degree | foreign.born | LTHS +
                   some.college + college + advanced.degree"), 
                 formLM = formula("age ~ ln.real.wage + female + 
                   some.college + college + advanced.degree + foreign.born + LTHS +
                   some.college + college + advanced.degree"), 
                 formGAM = formula("age ~ s(ln.real.wage) + female + 
                   some.college + college + advanced.degree + foreign.born + LTHS +
                   some.college + college + advanced.degree"), 
                 variables = c("age", "female", "college"),
                 kind = "barplot",
                 heapratio = 0.27, 
                 seed = 911){
  
  if(!is.null(seed)) set.seed(seed)
  library(data.table)
  df_heap <- df
  df_heap <- introduce_heaping(df, "age", heapratio, interval = 5)
  df_after_model <- df_after <- df_heap
  df_after$age <- correctHeaps2(df_heap$age, heaps = "5year", method = "lnorm")
  df_after_model$age <- correctHeaps2(df_heap$age, heaps = "5year", method = "lnorm",
                                           model = formLM, dataModel = df_heap)
  SEQ <- seq(min(df_heap$age, na.rm = TRUE), max(df_heap$age, na.rm = TRUE), 1)
  hist_data <- hist(df$age, plot = FALSE, breaks = SEQ)  # Simulate some histogram data
  # Extract bin centers and counts
  bin_centers <- hist_data$mids
  counts <- hist_data$counts
  # # Fit a GAM model to smooth the histogram data
  # gam_model <- gam(counts ~ s(bin_centers))
  # # Predict smoothed counts
  # df_gam <- round(predict(gam_model, newdata = data.frame(bin_centers = SEQ)))
  # 
  
  df_tmp <- data.frame(x = bin_centers, y = counts)
  mod_gam <- gam(y ~ s(x), data = df_tmp)
  df_gam <- round(predict(mod_gam, newdata = df_tmp, type = "response"))
  names(df_gam) <- names(table(df$age))
  
  
  res <- NULL
  
  if(kind == "barplot"){
    t1 <- table(df$age)
    t2 <- table(df_heap$age)
    t3 <- table(df_after$age)
    t4 <- table(df_after_model$age)
    t5 <- df_gam
    xlims <- max(c(t1, t2, t3, t4, t5)) 
    par(mfrow = c(2,3))
    barplot(t1, ylim = c(0, xlims), xlab = "age (original)")
    barplot(t2, ylim = c(0, xlims), xlab = "age (heaping)")
    barplot(t3, ylim = c(0, xlims), xlab = "age (corrected)")
    barplot(t4, ylim = c(0, xlims), xlab = "age (corrected + model)")
    barplot(t5, ylim = c(0, xlims), xlab = "age (smoothed gam)")
  }
  
  if(kind == "densityplot"){
    d <- density(df$age, bw = 1)
    par(mar = c(4,4,0.3,0.1))
    plot(d, main = "", xlab = "age", lwd = 2)
    lines(density(df_heap$age, bw = d$bw), col = "red", lwd = 2)
    lines(density(df_after$age, bw = d$bw), col = "blue")
    lines(density(df_after_model$age, bw = d$bw), col = "magenta")
    browser()
    lines(y = df_gam / sum(df_gam), x = as.numeric(names(df_gam)), col = "grey")
    legend("topright", legend = c("original", "heaping", "corrected", "corrected (model)", "smoothed (gam)"),
           col = c("black", "red", "blue", "magenta", "grey"), lwd = c(2,2,1,1,1))
  }
  
  if(kind %in% c("twofold", "beta.diff")){
    results <- oaxaca(formula = form, data = df, R = 1000)
    results_heap <- oaxaca(formula = form, data = df_heap, R = 1000)
    results_after <- oaxaca(formula = form, data = df_after, R = 1000)
    results_after_model <- oaxaca(formula = form, data = df_after_model, R = 1000)
  }
  
  if(kind == "twofold"){
    columns <- c("group.weight", "coef(unexplained A)", "coef(unexplained B)")
    res <- list("df" = results$twofold$variables[[5]][variables, columns],
      "df_heap" = results_heap$twofold$variables[[5]][variables, columns],
      "df_after" = results_after$twofold$variables[[5]][variables, columns],
      "df_after_model" = results_after_model$twofold$variables[[5]][variables, columns])
    res <- lapply(res, data.frame)
    res <- data.table::rbindlist(res,  use.names=TRUE, fill=TRUE, idcol = "ID")
    res <- res[seq(1, 12, 3), -2]
    colnames(res) <- c("data", "coef unexplained A", "coef unexplained B")
    res$data <- c("original", "heaping", "after", "after with model")
    res$diffA <- res[, 2]
    res$diffB <- res[, 3]
    res$diffA <- c(0, 100 * (res$`coef unexplained A`[1] - res$`coef unexplained A`[2]) / res$`coef unexplained A`[1], 
                  100 * (res$`coef unexplained A`[1] - res$`coef unexplained A`[3]) / res$`coef unexplained A`[1],
                  100 * (res$`coef unexplained A`[1] - res$`coef unexplained A`[4]) / res$`coef unexplained A`[1])
    res$diffB <- c(0, 100 * (res$`coef unexplained B`[1] - res$`coef unexplained B`[2]) / res$`coef unexplained B`[1], 
                  100 * (res$`coef unexplained B`[1] - res$`coef unexplained B`[3]) / res$`coef unexplained B`[1],
                  100 * (res$`coef unexplained B`[1] - res$`coef unexplained B`[4]) / res$`coef unexplained B`[1])
    res[, 2:ncol(res)] <- round(res[, 2:ncol(res)], 3)
  }
  
  if(kind == "beta.diff"){
   res <- c(results$beta$beta.diff["age"],
    results_heap$beta$beta.diff["age"],
    results_after$beta$beta.diff["age"],
    results_after_model$beta$beta.diff["age"])
   names(res) <- c("original", "heaping", "after", "after with model")
  }
  
  
  
  if(kind %in% c("coef", "conf.int")){
    lm1 <- lm(formula = formLM, data = df)
#    coef(lm1)[2]
    
    lm1_heap <- lm(formula = formLM, data = df_heap)

    lm1_after <- lm(formula = formLM, data = df_after)

    lm1_after_model <- lm(formula = formLM, data = df_after_model)
  }
  
  if(kind == "coef"){
    res <- c(coef(lm1)[2], coef(lm1_heap)[2], coef(lm1_after)[2], coef(lm1_after_model)[2])
    names(res) <- c("original", "heaping", "after", "after with model")
  }
  
  if(kind == "conf.int"){
    res <- list("df" = confint(lm1)[2, ],
                          "df_heap" = confint(lm1_heap)[2, ],
                          "df_after" = confint(lm1_after)[2, ],
                          "df_after_model" = confint(lm1_after_model)[2, ])
    #res$kind <- c("original", "heaping", "after", "after with model")
    res <- lapply(res, matrix, ncol = 2)
    res <- lapply(res, data.frame)
    res <- data.frame(rbindlist(res))
    rownames(res) <- c("original", "heaping", "after", "after with model")
  }

  if(kind %in% c("RMSE", "Rsquared", "MAE", "RMSESD", "MAESD", "all_caret")){
    require(caret)
    df$age <- as.numeric(as.character(cut(df$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    df_heap$age <- as.numeric(as.character(cut(df_heap$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    df_after$age <- as.numeric(as.character(cut(df_after$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    df_after_model$age <- as.numeric(as.character(cut(df_after_model$age, breaks=c(1,9,19,29,39,49,59,69,100), labels=1:8)))
    # Define the formula for the model
    # formula <- agecut ~ ln.real.wage + 
    #   female + 
    #   LTHS * foreign.born * advanced.degree + 
    #   some.college + 
    #   college
    # Define the training control with 10-fold cross-validation and 3 repetitions
    train_control <- trainControl(method = "repeatedcv", 
                                  number = 10, 
                                  repeats = 3)
    lm_df_age_caret <- train(formLM, 
                                  data = na.omit(df), 
                                  method = "lm", 
                                  trControl = train_control)
    lm_df_heap_age_caret <- train(formLM, 
                                       data = na.omit(df_heap), 
                                       method = "lm", 
                                       trControl = train_control)
    lm_df_after_age_caret <- train(formLM, 
                                        data = na.omit(df_after), 
                                        method = "lm", 
                                        trControl = train_control)
    lm_df_after_model_age_caret <- train(formLM, 
                                              data = na.omit(df_after_model), 
                                              method = "lm", 
                                              trControl = train_control)

  }
  if(kind == "RMSE"){
    res <- c(lm_df_age_caret$results$RMSE, 
             lm_df_heap_age_caret$results$RMSE, 
             lm_df_after_age_caret$results$RMSE, 
             lm_df_after_model_age_caret$results$RMSE)
    names(res) <- c("original", "heaping", "after", "after with model")
  }
  if(kind == "Rsquared"){
    res <- c(lm_df_age_caret$results$Rsquared, 
             lm_df_heap_age_caret$results$Rsquared, 
             lm_df_after_age_caret$results$Rsquared, 
             lm_df_after_model_age_caret$results$Rsquared)
    names(res) <- c("original", "heaping", "after", "after with model")
  }
  if(kind == "MAE"){
    res <- c(lm_df_age_caret$results$MAE, 
             lm_df_heap_age_caret$results$MAE, 
             lm_df_after_age_caret$results$MAE, 
             lm_df_after_model_age_caret$results$MAE)
    names(res) <- c("original", "heaping", "after", "after with model")
  }
  if(kind == "RMSESD"){
    res <- c(lm_df_age_caret$results$RMSESD, 
             lm_df_heap_age_caret$results$RMSESD, 
             lm_df_after_age_caret$results$RMSESD, 
             lm_df_after_model_age_caret$results$RMSESD)
    names(res) <- c("original", "heaping", "after", "after with model")
  }
  if(kind == "MAESD"){
    res <- c(lm_df_age_caret$results$MAESD, 
             lm_df_heap_age_caret$results$MAESD, 
             lm_df_after_age_caret$results$MAESD, 
             lm_df_after_model_age_caret$results$MAESD)
    names(res) <- c("original", "heaping", "after", "after with model")
  }
  if(kind == "all_caret"){
    res <- rbind(lm_df_age_caret$results[, -1], lm_df_heap_age_caret$results[, -1], 
          lm_df_after_age_caret$results[, -1],  lm_df_after_model_age_caret$results[, -1])
    rownames(res) <- c("original", "heaping", "after", "after with model")
    res <- as.matrix(res)
  }
  return(res)
}

pdf("figures/chicago_barplot.pdf")
chic(chicago)
dev.off()
par(mfrow = c(1,1))
pdf("figures/chicago_denplot.pdf")
chic(chicago, kind = "densityplot")
dev.off()


tf <- chic(chicago, kind = "twofold")
library(xtable)
xtable(tf, digits = 3, label = "tab:twofold",caption = "Twofold Blinder-Oaxaca decomposition. Coefficients from a regression on
observations from Groups A and B, respectively, in the reference coefficient vector, i.e. the set of
regression coefficients that would emerge in a world of no labor market discrimination.")

R <- 10

set.seed(11)
tf_rep <- replicate(R, as.matrix(chic(chicago, kind = "twofold", seed = NULL)[, 2:5]))
tf_rep


bd <- chic(chicago, kind = "beta.diff")
bd

set.seed(11)
bd_rep <- replicate(R, chic(chicago, kind = "beta.diff", seed = NULL))
bd_rep

c1 <- chic(chicago, kind = "coef")
c1
set.seed(11)
c1_rep <- replicate(R, chic(chicago, kind = "coef", seed = NULL))
c1_rep

c1_confint <- chic(chicago, kind = "conf.int")
c1_confint

set.seed(11)
c1_confint_rep <- replicate(R, as.matrix(chic(chicago, kind = "conf.int", seed = NULL)))
c1_confint_rep

cc <- chic(chicago, kind = "all_caret")
set.seed(11)
cc_rep <- replicate(R, chic(chicago, kind = "all_caret", seed = NULL))
cc_rep


library(dplyr)
data(eusilc13puf, package = "simPop")
eusilc13puf$age <- as.numeric(as.character(eusilc13puf$age)) + 1
eusilc13puf$linc <- log(eusilc13puf$pgrossIncome + 1)
eusilc13puf$citizenship <- ifelse(eusilc13puf$pb220a == "AT", TRUE, FALSE)
form1 <- formula("linc ~ age +
                rb090 + 
                pl031 + # Derzeitige Hauptaktivität
                pb190 + # Familienstand + 
                pe040 + 
                pl111 | citizenship")
form2 <-  formula("age ~ linc + 
                rb090 + 
                pl031 + # Derzeitige Hauptaktivität
                pb220a + 
                pb190 + # Familienstand + 
                pe040 + # Höchster Bildungsabschluss 
                pl111")
eusilc13puf <- handle_missing_values(eusilc13puf)

pdf("figures/eusilc_barplot.pdf")
chic(eusilc13puf, form = form1, formLM = form2, seed = 11, heapratio = 0.15)
dev.off()
pdf("figures/eusilc_densityplot.pdf")
chic(eusilc13puf, form = form1, formLM = form2, kind = "densityplot")
dev.off()
#chic(eusilc13puf, form = form1, formLM = form2, kind = "twofold")
#chic(eusilc13puf, form = form1, formLM = form2, kind = "beta.diff")
silc_coef <- chic(eusilc13puf, form = form1, formLM = form2, kind = "coef")
set.seed(11)
silc_coef_rep <- replicate(R, chic(eusilc13puf, form = form1, formLM = form2, kind = "coef", seed = NULL))

silc_confint <- chic(eusilc13puf, form = form1, formLM = form2, kind = "conf.int")
set.seed(11)
silc_confint_rep <- replicate(R, chic(eusilc13puf, form = form1, formLM = form2, kind = "conf.int", seed = NULL))

silc_cc <- chic(eusilc13puf, form = form1, formLM = form2, kind = "all_caret")
set.seed(11)
silc_cc_rep <- replicate(R, chic(eusilc13puf, form = form1, formLM = form2, kind = "all_caret", seed = NULL))


##########################################################
### old ##################################################


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
