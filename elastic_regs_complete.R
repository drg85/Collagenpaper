library(metafor)
library(readxl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
datain <- read_excel("allelasticfinal.xlsx")


# data in
dat <- data.frame(
  study = datain$study,
  n1 = datain$NE,
  mean1 = datain$ME,
  sd1 = datain$SE,
  n2 = datain$NC,
  mean2 = datain$MC,
  sd2 = datain$SC,
  dose = datain$dose,       # grams/day
  weeks = datain$Time,
  pressure = datain$Pressure,
  suction = datain$Suction,
  ma1 = datain$ma1,
  ma2 = datain$ma2
)
dat$total_dose <- dat$dose*dat$weeks*7


#Set to 1 for just suction methods, 2 for just pressure
methods <- 0
if(methods ==1){
  dat <- dat[dat$suction==1,]
}

if(methods ==2){
  dat <- dat[dat$pressure==1,]
}

#if you just want De Miranda, set this to 1
mir <- 0
if(mir==1){
  dat <- dat[dat$ma1 ==1,]
}

# Compute raw mean differences and variances
escalc_data <- escalc(
  measure = "SMD",
  m1i = mean1, sd1i = sd1, n1i = n1,
  m2i = mean2, sd2i = sd2, n2i = n2,
  data = dat
)

# Fit model with total dose and duration (weeks)
res_weekly <- rma(yi, vi, mods = ~ total_dose + weeks , data = escalc_data)


#Daily variant 
res_daily <- rma(yi, vi, mods = ~ dose + weeks , data = escalc_data)

#Dose variant 
res_dose <- rma(yi, vi, mods = ~ dose , data = escalc_data)
# Predict values across a sequence of daily doses (ignoring weeks)
new_dose <- seq(min(dat$dose), max(dat$dose), length.out = 100)
preds <- predict(res_dose, newmods = new_dose)  # res3 uses mods = ~ dose + weeks

bubble_sizes <- sqrt(1 / escalc_data$vi)
bubble_sizes <- pmax(bubble_sizes, 1.5)  # Set minimum bubble size to 1.5

# Plot raw effect sizes vs dose
plot(dat$dose, escalc_data$yi,
     xlab = "Daily Dose (g/day)",
     ylab = "Mean Difference",
     pch = 21, bg = "lightblue",
     cex =  bubble_sizes,  # Bubble size ~ precision
     main = "Elastic meta-Regression"
)

# Add regression line
lines(new_dose, preds$pred, col = "red", lwd = 2)

# Add confidence interval bands
lines(new_dose, preds$ci.lb, col = "red", lty = 2)
lines(new_dose, preds$ci.ub, col = "red", lty = 2)



library(metafor)

# Fit meta-regression with cumulative dose
res_cum <- rma(yi, vi, mods = ~ total_dose, data = escalc_data)

# Predict across a sequence of cumulative doses
new_total <- seq(min(dat$total_dose), max(dat$total_dose), length.out = 100)
preds_cum <- predict(res_cum, newmods = new_total)

# Bubble sizes ~ precision
bubble_sizes <- sqrt(1 / escalc_data$vi)
bubble_sizes <- pmin(bubble_sizes, 5)  # Limit max size

# Plot
plot(dat$total_dose, escalc_data$yi,
     xlab = "Cumulative Dose (g)",
     ylab = "Standardized Mean Difference (SMD)",
     pch = 21, bg = "lightblue",
     cex = bubble_sizes,
     main = "(b) Elasticity Meta-Regression vs Cumulative Dose"
)

# Add shaded confidence band
polygon(c(new_total, rev(new_total)),
        c(preds_cum$ci.lb, rev(preds_cum$ci.ub)),
        col = rgb(1, 0, 0, 0.2), border = NA)

# Add regression line
lines(new_total, preds_cum$pred, col = "red", lwd = 2)

