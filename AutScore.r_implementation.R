# Implementing AutScore.r...
start_time <- Sys.time()  # Record start time

library(tidyverse)
library(readr)
library(forcats)

# Import the list of LP/P/LGDs
data = read_csv("C:/Users/au759264/OneDrive - Aarhus universitet/Postdoc works/PhD_Codes/AutScore_paper_submission files-Scientific reports/2nd Revision/github_data_upload.csv")

# Preprocessing
data2=data %>%
  mutate(
    AutScore_Cat=ifelse(AutScore > 11, 1, 0), 
    AutoCasC_Cat=ifelse(AutoCasC > 6, 1, 0), 
    D_Aut_CasC=AutScore_Cat-AutoCasC_Cat, 
    D_CasC_Noa=AutoCasC_Cat-Clinical_Score, 
    D_Aut_Noa=AutScore_Cat-Clinical_Score,
    SFARI = as.factor(ifelse(SFARI > 1, 2, SFARI)),
    SFARI = fct_relevel(SFARI, "0"),
    GDA_Score = as.numeric(GDA_Score),
    GDA_C=as.factor(ifelse(GDA_Score < 0.10, 0,
                           ifelse(GDA_Score >= 0.10 & GDA_Score < 0.30, 1,
                                  ifelse(GDA_Score >= 0.30 & GDA_Score < 0.49, 2, 3)))),
    GDA_C = fct_relevel (GDA_C, "0"),
    I_C=as.factor(ifelse(I == 6, "P",
                         ifelse(I == 3, "LP", "VUS"))),
    I_C = fct_relevel (I_C, "VUS"),
    D_C=as.factor(ifelse(D==1, 2,
                         ifelse(D==2, 2, 
                                ifelse(D==-1, 0,
                                       ifelse(D==-2, 0, 1))))),
    D_C = fct_relevel (D_C, "0"),
    C_C=as.factor(ifelse(C==1, 1,
                         ifelse(C==3, 1, 0))),
    C_C = fct_relevel(C_C, "0"),
    H_C=ifelse(H > 0, 1, 0),
    P = as.factor(ifelse(P <= 3, 0,
                         ifelse(P > 3 & P <= 5, 1, 2))),
    P = fct_relevel (P, "0")
  ) %>%
  select(CHROM, REF, ALT, POS = `POS-hg38`, Consequences = Conseq_Revised, Gene_Type, 
         Clinical_Score, GDA_C, I_C, Psi = P, D_C, SFARI, H, H_C, C, C_C, 
         AutScore, AutoCasC)

# Model fitting
model = glm(Clinical_Score ~ I_C + GDA_C + Psi + D_C + SFARI + H_C + C_C, 
            data = data2, family = binomial)

# Compute predicted probabilities
predicted_probs = predict(model, type = "response")

# Store predicted probabilities in the data frame
data2 = data2 %>%
  mutate(AutScore.r = round(predicted_probs, 3)) %>%
  select(CHROM, REF, ALT, POS, Consequences, Gene_Type, 
         Clinical_Score, AutScore, AutoCasC, AutScore.r)

# Export the AutScore.r output
write.csv(data2, "C:/Users/au759264/OneDrive - Aarhus universitet/Postdoc works/PhD_Codes/AutScore_paper_submission files-Scientific reports/2nd Revision/AutScore.r.output.csv")

# Extract the model summary
summary_model = summary(model)

# Extract the coefficients (log-odds) from the summary
coef = summary_model$coefficients

# Extract the coefficients (log-odds)
round((coef(summary(model))), 3)

# Compute Odds Ratios (OR) by exponentiation the coefficients
round((exp(coef[, "Estimate"])), 2)

# Compute 95% CI for the OR
round((exp(coef[, "Estimate"] - 1.96 * coef[, "Std. Error"])), 2)
round((exp(coef[, "Estimate"] + 1.96 * coef[, "Std. Error"])), 2)


# Load the package
library(pROC)

# Create the ROC curve
roc_curve = roc(data2$Clinical_Score, data2$AutScore.r)

# Plot the ROC curve

plot(roc_curve, main = NA, col = "blue", legacy.axes = TRUE)

# Calculate AUC
auc_value = auc(roc_curve)
print(paste("AUC:", auc_value))

# Calculate the 95% confidence interval for AUC
ci_auc = ci.auc(roc_curve)
print(paste("95% CI for AUC:", ci_auc))

# Find the optimal threshold using Youden's index
optimal_cutoff = coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal cut-off point
print(optimal_cutoff)


# Get sensitivity and specificity at the optimal cut-off
sens_spec = coords(roc_curve, "best", ret = c("sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(sens_spec)


end_time = Sys.time()  # Record end time

# Compute and print runtime
runtime = end_time - start_time
print(paste("Total runtime:", runtime)) # 2.178 seconds
