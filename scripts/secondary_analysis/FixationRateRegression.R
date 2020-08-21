

### Project: Poliovirus sequencing from Matlab
### Purpose: Beta regression of gatekeeper mutations
### Working directory: Poliovirus_Intrahost

# ================== Import packages, read in data =================

library(tidyverse)
library(wesanderson)
library(ggbeeswarm)
library(reshape2)
library(lubridate)
library(cowplot)
library(patchwork)
library(betareg)

frequencies <- read.csv("data/processed/gatekeepers_cross_section.csv")

palette <- wes_palette("Darjeeling1")
darjeeling2_palette <- wes_palette("Darjeeling2")

# ================================ Support functions ================================

# Cross-sectional plot
PlotCrossSection <- function(data, data_box, mut, colors)
{
  
  plot <- ggplot() + 
    geom_boxplot(data = filter(data_box, mutation == mut), aes(x = as.factor(collectionWeekRelative), y = frequency), notch = FALSE, outlier.shape = NA, alpha = 0.9) + 
    theme_bw() + 
    geom_jitter(data = filter(data, mutation == mut), aes(x = as.factor(collectionWeekRelative), y = frequency, fill = mutation), width = 0.2) +
    xlab("Week Post Vaccination") + 
    ylab("Frequency") + 
    scale_fill_manual(name = "Mutation", values = colors) +
    theme(legend.position = "none")
  
  return(plot)
}

# ================================ Beta regression =========================

frequencies_beta <- mutate(frequencies, frequency_beta = ifelse(frequency == 0, 0.0000001, frequency))
frequencies_beta <- mutate(frequencies_beta, frequency_beta = ifelse(frequency_beta == 1, 0.9999999, frequency_beta))

betareg.model.481 <- betareg(frequency_beta ~ collectionWeekRelative, data = filter(frequencies_beta, mutation == "A481G" & collectionWeekRelative <= 4), link = "logit")
betareg.model.2909 <- betareg(frequency_beta ~ collectionWeekRelative, data = filter(frequencies_beta, mutation == "U2909C" & collectionWeekRelative <= 4), link = "logit")
betareg.model.398 <- betareg(frequency_beta ~ collectionWeekRelative, data = filter(frequencies_beta, mutation == "U398C" & !collectionWeekRelative %in% c(1, 2)), link = "logit")

# Plot model predictions further in time.
weeks <- seq(0, 100, by = 1)
prediction_df_full <- data.frame(collectionWeekRelative = weeks)
prediction_df_full$predicted_481 <- predict(betareg.model.481, newdata = prediction_df_full) # Week where it exceeds 0.95: week 6
prediction_df_full$predicted_2909 <- predict(betareg.model.2909, newdata = prediction_df_full) # Week where it exceeds 0.95: week 13
prediction_df_full$predicted_398 <- predict(betareg.model.398, newdata = prediction_df_full) # Week where it exceeds 0.95: week 46

beta.predict.plot <- ggplot() +
  geom_line(data = prediction_df_full, aes(x = collectionWeekRelative, y = predicted_481)) +
  geom_line(data = prediction_df_full, aes(x = collectionWeekRelative, y = predicted_2909)) +
  geom_line(data = prediction_df_full, aes(x = collectionWeekRelative, y = predicted_398)) +
  ylab("Predicted Frequency") +
  xlab("Week Post Vaccination") +
  xlim(0, 50)

# Get alpha and beta from each model.
estBetaParams <- function(mu, var) 
{
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams(as.numeric(coef(betareg.model.481)[2]), as.numeric(coef(betareg.model.481)[3]))
estBetaParams(as.numeric(coef(betareg.model.2909)[2]), as.numeric(coef(betareg.model.2909)[3]))
estBetaParams(as.numeric(coef(betareg.model.398)[2]), as.numeric(coef(betareg.model.398)[3]))

# ================================= Plot with data ================================

weeks <- seq(1, 6, by = 1)
prediction_df <- data.frame(collectionWeekRelative = weeks)
prediction_df$predicted_481 <- predict(betareg.model.481, newdata = prediction_df)
prediction_df$predicted_2909 <- predict(betareg.model.2909, newdata = prediction_df)
prediction_df$predicted_398 <- predict(betareg.model.398, newdata = prediction_df)

frequencies %>%
  group_by(mutation, collectionWeekRelative) %>%
  filter(n() > 4) -> complete_data_all_box

cross.sectional.481.plot <- PlotCrossSection(frequencies, complete_data_all_box, "A481G", c(palette[4])) + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481, color = "red")) + 
  ggtitle("A481G")
cross.sectional.2909.plot <- PlotCrossSection(frequencies, complete_data_all_box, "U2909C", c(darjeeling2_palette[2])) + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909, color = "red")) +
  ggtitle("U2909C")
cross.sectional.398.plot <- PlotCrossSection(frequencies, complete_data_all_box, "U398C", c(palette[5])) + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398, color = "red")) +
  ggtitle("U398C")

cross.section.plot <- plot_grid(cross.sectional.481.plot, cross.sectional.2909.plot, cross.sectional.398.plot, ncol = 3) # 10 by 3.5

# ======================================== Exponential nls fitting ====================================

model.481 <- nls(frequency ~ 1 - exp(-k*collectionWeekRelative), data = filter(frequencies, mutation == "A481G"), start = list(k = 1), trace = TRUE)
model.2909 <- nls(frequency ~ 1 - exp(-k*collectionWeekRelative), data = filter(frequencies, mutation == "U2909C"), start = list(k = 1), trace = TRUE)
model.398 <- nls(frequency ~ 1 - exp(-k*collectionWeekRelative), data = filter(frequencies, mutation == "U398C"), start = list(k = 1), trace = TRUE)

# Plot the predictions with data.
weeks <- seq(1, 6, by = 1)
prediction_df <- data.frame(collectionWeekRelative = weeks)
prediction_df$predicted_481 <- predict(model.481, newdata = prediction_df)
prediction_df$predicted_2909 <- predict(model.2909, newdata = prediction_df)
prediction_df$predicted_398 <- predict(model.398, newdata = prediction_df)
prediction_df <- mutate(prediction_df, predicted_481_IDM = 1 - exp(-(1/7)*collectionWeekRelative))
prediction_df <- mutate(prediction_df, predicted_2909_IDM = 1 - exp(-(1/21)*collectionWeekRelative))
prediction_df <- mutate(prediction_df, predicted_398_IDM = 1 - exp(-(1/47)*collectionWeekRelative))
prediction_df$predicted_481_beta <- predict(betareg.model.481, newdata = prediction_df)
prediction_df$predicted_2909_beta <- predict(betareg.model.2909, newdata = prediction_df)
prediction_df$predicted_398_beta <- predict(betareg.model.398, newdata = prediction_df)

cross.sectional.481.nls.plot <- PlotCrossSection(frequencies, complete_data_all_box, "A481G", c(palette[4])) + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481), color = "red") + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481_IDM), color = "blue") + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481_beta), color = "purple") + 
  ggtitle("A481G")
cross.sectional.2909.nls.plot <- PlotCrossSection(frequencies, complete_data_all_box, "U2909C", c(darjeeling2_palette[2])) + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909), color = "red") +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909_IDM), color = "blue") + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909_beta), color = "purple") + 
  ggtitle("U2909C")
cross.sectional.398.nls.plot <- PlotCrossSection(frequencies, complete_data_all_box, "U398C", c(palette[5])) + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398, color = "red")) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398_IDM), color = "blue") + 
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398_beta), color = "purple") + 
  ggtitle("U398C")

cross.section.plot.nls <- plot_grid(cross.sectional.481.nls.plot, cross.sectional.2909.nls.plot, cross.sectional.398.nls.plot, ncol = 3) # 10 by 3.5

# Plot model predictions further in time.
weeks <- seq(0, 100, by = 1)
prediction_df <- data.frame(collectionWeekRelative = weeks)
prediction_df$predicted_481 <- predict(model.481, newdata = prediction_df)
prediction_df$predicted_2909 <- predict(model.2909, newdata = prediction_df)
prediction_df$predicted_398 <- predict(model.398, newdata = prediction_df)

prediction_df <- mutate(prediction_df, predicted_481_IDM = 1 - exp(-(1/7)*collectionWeekRelative*7))
prediction_df <- mutate(prediction_df, predicted_2909_IDM = 1 - exp(-(1/21)*collectionWeekRelative*7))
prediction_df <- mutate(prediction_df, predicted_398_IDM = 1 - exp(-(1/47)*collectionWeekRelative*7))

prediction_df$predicted_481_beta <- predict(betareg.model.481, newdata = prediction_df)
prediction_df$predicted_2909_beta <- predict(betareg.model.2909, newdata = prediction_df)
prediction_df$predicted_398_beta <- predict(betareg.model.398, newdata = prediction_df)

a = 1
nls.predict.plot.481 <- ggplot() +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481), color = "red", alpha = a) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481_IDM), color = "blue", alpha = a) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_481_beta), color = "purple", alpha = a) + 
  ylab("Predicted Frequency") +
  xlab("Week Post Vaccination") +
  xlim(c(0, 50))

nls.predict.plot.2909 <- ggplot() +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909), color = "red", alpha = a) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909_IDM), color = "blue", alpha = a) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_2909_beta), color = "purple", alpha = a) + 
  ylab("Predicted Frequency") +
  xlab("Week Post Vaccination") +
  xlim(c(0, 50))

nls.predict.plot.398 <- ggplot() +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398), color = "red", alpha = a) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398_IDM), color = "blue", alpha = a) +
  geom_line(data = prediction_df, aes(x = collectionWeekRelative, y = predicted_398_beta), color = "purple", alpha = a) + 
  ylab("Predicted Frequency") +
  xlab("Week Post Vaccination") +
  xlim(c(0, 50))

nls.predict.plot.all <- plot_grid(nls.predict.plot.481, nls.predict.plot.2909, nls.predict.plot.398, ncol = 3)
