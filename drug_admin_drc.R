library(tidyverse)
library(drc)
library(ggthemr)

df <- read.csv('1.29.26.csv')

###############################################################################
#                         ***  HELPER FUNCTIONS  ***
###############################################################################

# Set ggplot2 theme
ggthemr('earth')

# Add a norm_counts column at the end of the original .csv file

low_dose_tbl <- df |> 
  filter(dose == 0.00001) |> 
  group_by(cell_line) |> 
  summarise(
    baseline = mean(counts, na.rm = T),
    .groups = 'drop'
  )

df_norm <- df |> 
  left_join(low_dose_tbl, by = "cell_line") |> 
  mutate(
    norm_counts = counts / baseline
  )

df <- subset(df_norm, select=-baseline)
rm(list = c('df_norm', 'low_dose_tbl'))

my_palette <- c(
  "55-3"    = "blue",
  "55-7"    = "red",
  "57-2"    = "green",
  "par" = "purple"
)

# 1. Create 100 points evenly spaced in log10‑space to help smooth the line later
log_grid <- seq(
  from = log10(min(df$dose)),   # log10(0.0001) = -4
  to   = log10(max(df$dose)),   # log10(1000)   =  3
  length.out = 100
)

# 2. Transform back to linear scale
grid_linear <- 10^log_grid

# Create a data frame of doses spanning the observed range (log-scale for even spacing)
new_doses <- data.frame(Dose = grid_linear)

# FUNCTION TO CREATE LL.4 DOSE RESPONSE MODEL

doseResponseModelViability <- function(df, cellLine) {
  cellLine_df <- df %>% filter(cell_line == cellLine)
  model_cellLine <- drm(
    formula = viability ~ dose,
    data = cellLine_df,
    fct = LL.4()
  )
  
  return(model_cellLine)
}


doseResponseModelCC <- function(df, cellLine) {
  cellLine_df <- df %>% filter(cell_line == cellLine)
  model_cellLine <- drm(
    formula = norm_counts ~ dose,
    data = cellLine_df,
    fct = LL.4()
  )
  
  return(model_cellLine)
}

# Create Prediction Response for eventual line-of-best fit

predictedResponse <- function(model, new_doses) {
  predicted_cellLine <- data.frame(
    Dose = new_doses$Dose,
    Response = predict(model, newdata = new_doses)
  )
}

###############################################################################
#                     ***  ANALYSIS WITH VIABILITY  ***
###############################################################################

model_553 <- doseResponseModelViability(df, "55-3")
predicted_553 <- predictedResponse(model_553, new_doses)
# ED(model, 50) = dose for 50% response; interval = "fls" for confidence limits
ec50_553 <- ED(model_553, 50, interval = "fls")  

# Predict responses for these doses using the model

model_557 <- doseResponseModelViability(df, "55-7")
predicted_557 <- predictedResponse(model_557, new_doses)
ec50_557 <- ED(model_557, 50, interval = "fls")  

model_572 <- doseResponseModelViability(df, "57-2")
predicted_572 <- predictedResponse(model_572, new_doses)
ec50_572 <- ED(model_572, 50, interval = "fls")  

model_par <- doseResponseModelViability(df, "par")
predicted_par <- predictedResponse(model_par, new_doses)
ec50_par <- ED(model_par, 50, interval = "fls")  

###############################################################################
#                         *** PLOTTING VIABILITY***
###############################################################################

pred_via_long <- bind_rows(
  predicted_553 %>% mutate(cell_line = "55-3"),
  predicted_557 %>% mutate(cell_line = "55-7"),
  predicted_572 %>% mutate(cell_line = "57-2"),
  predicted_par %>% mutate(cell_line = "par")
)

ggplot() +
  
  geom_point(
    data = df,
    aes(x = dose, y = viability, colour = cell_line),
    size = 2
  ) +
  
  ## 4️⃣ Predicted curves – also coloured by cell_line
  geom_line(
    data = pred_via_long,
    aes(x = Dose, y = Response, colour = cell_line),
    linewidth = 1
  ) +
  
  ## 5️⃣ Keep the log‑scaled x‑axis (your custom breaks/labels)
  scale_x_log10(
    breaks = c(0.0001, 0.1, 1, 10, 100, 1000),
    labels = c("0.0001", "0.1", "1", "10", "100", "1000")
  ) +
  
  ## 6️⃣ Apply the manual palette (this creates the legend automatically)
  scale_colour_manual(
    name = "Cell line",   # legend title
    values = my_palette,
    # optional: nicer labels in the legend
    labels = c(
      "55‑3"    = "55‑3",
      "55‑7"    = "55‑7",
      "57‑2"    = "57‑2",
      "par" = "Parental"
    )
  ) +
  
  labs(
    x = "Dose (nM)",
    y = "Cell Viability (%)",
    title = "Dexamethasone Dose Response Curve for Nalm 6"
  ) 


###############################################################################
#                   ***  ANALYSIS USING CELL COUNTS  ***
###############################################################################

model_cc_553 <- doseResponseModelCC(df, "55-3")
predicted_cc_553 <- predictedResponse(model_cc_553, new_doses)

model_cc_557 <- doseResponseModelCC(df, "55-7")
predicted_cc_557 <- predictedResponse(model_cc_557, new_doses)

model_cc_572 <- doseResponseModelCC(df, "57-2")
predicted_cc_572 <- predictedResponse(model_cc_572, new_doses)

model_cc_par <- doseResponseModelCC(df, "par")
predicted_cc_par <- predictedResponse(model_cc_par, new_doses)


###############################################################################
#                         *** PLOTTING VIABILITY***
###############################################################################

# ------------------------------------------------------------------
# 1️⃣ Combine the four prediction data frames
# ------------------------------------------------------------------
pred_long <- bind_rows(
  predicted_cc_553 %>% mutate(cell_line = "55-3"),
  predicted_cc_557 %>% mutate(cell_line = "55-7"),
  predicted_cc_572 %>% mutate(cell_line = "57-2"),
  predicted_cc_par %>% mutate(cell_line = "par")
)

ggplot() +

  geom_point(
    data = df,
    aes(x = dose, y = norm_counts, colour = cell_line),
    size = 2
  ) +
  
  ## 4️⃣ Predicted curves – also coloured by cell_line
  geom_line(
    data = pred_long,
    aes(x = Dose, y = Response, colour = cell_line),
    linewidth = 1
  ) +
  
  ## 5️⃣ Keep the log‑scaled x‑axis (your custom breaks/labels)
  scale_x_log10(
    breaks = c(0.0001, 0.1, 1, 10, 100, 1000),
    labels = c("0.0001", "0.1", "1", "10", "100", "1000")
  ) +
  
  ## 6️⃣ Apply the manual palette (this creates the legend automatically)
  scale_colour_manual(
    name = "Cell line",   # legend title
    values = my_palette,
    # optional: nicer labels in the legend
    labels = c(
      "55‑3"    = "55‑3",
      "55‑7"    = "55‑7",
      "57‑2"    = "57‑2",
      "par" = "Parental"
    )
  ) +
  
  labs(
    x = "Dose (nM)",
    y = "Normalized Cell Counts",
    title = "Dexamethasone Dose Response Curve for Nalm 6"
  )
