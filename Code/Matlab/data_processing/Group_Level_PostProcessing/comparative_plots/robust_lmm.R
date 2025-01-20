library(R.matlab)

# Load the .mat file
print("Loading data_for_R.mat...")
mat_data <- readMat("data_for_R.mat") # Load entire file

# Debug: Check available variables
print(names(mat_data)) # Ensure 'data_struct' is listed

# Access the data_struct
if (!"data_struct" %in% names(mat_data)) {
  stop("data_struct not found in data_for_R.mat. Check MATLAB export.")
}
data_struct <- mat_data$data_struct

# Convert structure to a data frame
data <- as.data.frame(data_struct)

# Debug: Inspect data frame
print("Head of the data frame:")
print(head(data))
str(data)

# Debug: Inspect data frame
print("Head of the data frame:")
print(head(data))

# Ensure correct data types
data$Condition_ID <- as.factor(data$Condition_ID)
data$Subject_ID <- as.factor(data$Subject_ID)
data$IC_ID <- as.factor(data$IC_ID)

# Debug: Check factor levels
print("Levels of Condition_ID:")
print(levels(data$Condition_ID))

# Step 4: Fit the Robust LMM
library(robustlmm)
model <- rlmer(RMS_value ~ Condition_ID + (1 | Subject_ID) + (1 | Subject_ID:IC_ID), data = data)

# Debug: Print model summary
print("Model summary:")
print(summary(model))

# Step 5: Save Results
fixedEffects <- as.vector(fixef(model))
randomEffects <- as.vector(ranef(model)$Subject_ID)
robust_residuals <- residuals(model)

saveMat("results.mat", list(
  fixedEffects = fixedEffects,
  randomEffects = randomEffects,
  residuals = robust_residuals
))

print("Results saved to results.mat.")
