library(readr)
library(dplyr)

# Demo datasets for code development
if (FALSE) {
  prescriptions <- read_csv("Downloads/mimic-iv-clinical-database-demo-2.2/hosp/prescriptions.csv")
  patients <- read_csv("Downloads/mimic-iv-clinical-database-demo-2.2/hosp/patients.csv")
  omr <- read_csv("Downloads/mimic-iv-clinical-database-demo-2.2/hosp/omr.csv")
  d_labitems <- read_csv("Downloads/mimic-iv-clinical-database-demo-2.2/hosp/d_labitems.csv")
  labevents <- read_csv("Downloads/mimic-iv-clinical-database-demo-2.2/hosp/labevents.csv")
}

# preprocess the full datasets
if (FALSE) {
  # prescriptions
  prescriptions <- read_csv("Downloads/MIMIC/prescriptions.csv")
  Insulin <- prescriptions[ prescriptions[, "drug"] == "Insulin", c("subject_id", "starttime", "stoptime", "dose_val_rx")]
  write.csv(Insulin, "Downloads/MIMIC/insulin.csv")
  # glucose indicees
  lab_glucose_indices <- which(d_labitems$label == "Glucose") 
  # labevents
  labevents <- read_csv("Downloads/MIMIC/labevents.csv")
  Glucose <- labevents[labevents$itemid == "50931", c("subject_id", "charttime", "valuenum")]
  write.csv(Glucose, "Downloads/MIMIC/glucose.csv")
  # better in Bash: 
  # awk -F, 'NR==1 {print $2","$7","$10} $5=="50931" {print $2","$7","$10}'  Downloads/MIMIC/labevents.csv > Downloads/MIMIC/glucose.csv
}

Insulin <- read_csv("Downloads/MIMIC/insulin.csv")
patients <- read_csv("Downloads/MIMIC/patients.csv")
omr <- read_csv("Downloads/MIMIC/omr.csv")
d_labitems <- read_csv("Downloads/MIMIC/d_labitems.csv")
Glucose <- read_csv("Downloads/MIMIC/glucose.csv")

n <- nrow(Insulin)

Insulin$Age <- rep(NA, n)
Insulin$BMI <- rep(NA, n)
Insulin$Glucose_delta <- rep(NA, n)

for (i in 1:n) {
  patient_id <- as.numeric(Insulin[i, "subject_id"])
  
  index <- which(patients$subject_id == patient_id)
  Insulin$Age[i] <- patients$anchor_age[index]
  
  index <- sum(omr$subject_id == patient_id  & omr$result_name == "BMI (kg/m2)")
  if (index > 0){Insulin$BMI[i] <- omr$result_value[index]} # some NA and some fractions
  
  # find closest Glucose measurements
  Glucose_i <- Glucose[Glucose$subject_id == patient_id, ] 
  time_of_insulin <- Insulin$stoptime[i]
  
  time_before_i <- time_of_insulin - Glucose_i$charttime
  if (length(time_before_i[time_before_i > 0]) == 0) {next}
  index_time_before_i <- which(time_before_i == min(time_before_i[time_before_i > 0]))
  Glucose_before <- Glucose_i$valuenum[index_time_before_i]
  
  time_after_i <- Glucose_i$charttime - time_of_insulin
  if (length(time_after_i[time_after_i > 0]) == 0) {next}
  index_time_after_i <- which(time_after_i == min(time_after_i[time_after_i > 0]))
  Glucose_after <- Glucose_i$valuenum[index_time_after_i]
  
  if(Glucose_after > 0 & Glucose_before > 0 ) { Insulin$Glucose_delta[i] <- Glucose_after - Glucose_before}
}
# did not run completely through, but we got enough samples. # 1987 samples

Insulin <- Insulin[!is.na(Insulin$Glucose_delta), c("Age", "BMI", "Glucose_delta", "dose_val_rx")]
Insulin <- na.omit(Insulin)

# process the column BMI by computing the entries such as "106/63"
Insulin_clean <- Insulin %>%
  mutate(
    BMI = ifelse(
      grepl("/", BMI),                     
      sapply(strsplit(BMI, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      as.numeric(BMI)                      
    ),
    dose_val_rx = as.numeric(dose_val_rx)  
  )
Insulin_clean <- as.matrix(Insulin_clean)
Insulin_clean <- Insulin_clean[Insulin_clean[, 4] > 0, ]
colnames(Insulin_clean) <- NULL


# ------------------------------------------------------------------------------
# MIMIC Analysis
# ------------------------------------------------------------------------------


source("./Documents/papers/2025_sample_optimization/code/R/function_back_door.R")
data <- data.frame(C = rep(Inf, n), X_t = Insulin_clean[, 4], X_M = Insulin_clean[, 3])
data$X_B <- I(Insulin_clean[, 1:2])
data <- na.omit(data)
n <- nrow(data)

c1_fun <-  function(XBt) { rep(1, nrow(XBt))}
pi1_fun <- function(XBt) { rep(1, nrow(XBt)) }
indicator_2inf <- rep(TRUE, n)


set.seed(1234)
estimates <- estimate_parameters(data, pi1_fun)
mom <- compute_moments(data, estimates, pi1_fun)

# evaluate vanilla
avars_vanilla <- compute_oif_variance(data, estimates, pi1_fun, 
                                      pi1 = rep(1, dim(data)[1]))
cost_vanilla <- compute_expected_cost(rep(1, dim(data)[1]),
                                      c0 = 1, c1_fun(data[, c("X_B", "X_t")]),
                                      mom$Weight_2inf, indicator = indicator_2inf) # Why NA?
for (scale in seq(1.1, 2, by = 0.01)){
  print(scale)
  b0 <- cost_vanilla / scale
  opt <- compute_optimal_pi(data, estimates, 
                            c0 = 1, c1_fun = c1_fun,
                            b0 = b0, 
                            pi1_0_fun = pi1_fun, n_sub = 10000)
  avars_opt <- compute_oif_variance(data, estimates,
                                    pi1_fun,
                                    pi1 = opt$pi1_star)
  print(c(scale, avars_opt / avars_vanilla / scale))
  
}




