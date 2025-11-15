"""
Random Forest Epigenetic Clock - HPC Version
Predicts age from DNA methylation beta values using a Random Forest model
Includes SHAP analysis for feature importance
"""

# -----------------------------
# Load required packages
# -----------------------------
import pandas as pd
import numpy as np
import shap
import matplotlib
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, KFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import make_scorer, mean_absolute_error, mean_squared_error

# -----------------------------
# Load data
# -----------------------------
input_file = "subset_500_Beta.csv"
GSE55763_beta = pd.read_csv(input_file)
print("Dataset shape:", GSE55763_beta.shape)

# -----------------------------
# Prepare predictors and outcome
# -----------------------------
cpg_cols = [col for col in GSE55763_beta.columns if col.startswith("cg")]

X = GSE55763_beta[cpg_cols + ["gender"]].copy()
y = GSE55763_beta["age"]

X["gender"] = X["gender"].map({"M": 1, "F": 0})

# -----------------------------
# Split into training and testing sets
# -----------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=500
)

# Scale features
scaler = StandardScaler()
X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_test_scaled = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)

# -----------------------------
# Train Random Forest model
# -----------------------------
rf = RandomForestRegressor(n_estimators=300, random_state=500)
rf.fit(X_train_scaled, y_train)

# -----------------------------
# 5-Fold Cross-Validation
# -----------------------------
kf = KFold(n_splits=5, shuffle=True, random_state=500)
rf_cv = RandomForestRegressor(n_estimators=300, random_state=500)

scoring = {
    'R2': 'r2',
    'MAE': make_scorer(mean_absolute_error, greater_is_better=False),
    'RMSE': make_scorer(lambda y_true, y_pred: np.sqrt(mean_squared_error(y_true, y_pred)),
                         greater_is_better=False)
}

cv_results = cross_validate(
    rf_cv,
    X_train_scaled,
    y_train,
    cv=kf,
    scoring=scoring,
    return_train_score=False
)

print("Random Forest 5-Fold Cross-Validation Results:")
print(f"Mean R²:   {np.mean(cv_results['test_R2']):.2f}")
print(f"Mean MAE:  {abs(np.mean(cv_results['test_MAE'])):.2f}")
print(f"Mean RMSE: {abs(np.mean(cv_results['test_RMSE'])):.2f}")

# -----------------------------
# Metrics on Test Set
# -----------------------------
y_pred = rf.predict(X_test_scaled)
r2 = rf.score(X_test_scaled, y_test)
mae = mean_absolute_error(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))

print("\nTest Set Metrics:")
print(f"R²: {r2:.2f}")
print(f"MAE: {mae:.2f}")
print(f"RMSE: {rmse:.2f}")

# -----------------------------
# SHAP Analysis
# -----------------------------
explainer = shap.TreeExplainer(rf)

# Subset for faster SHAP computation (optional)
X_shap = X_test_scaled.sample(100, random_state=500)
shap_values = explainer.shap_values(X_shap)

plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_shap, plot_type="bar", color="coral", show=False)
plt.title("SHAP Summary Plot (Subset of Test Set)")
plt.tight_layout()
plt.savefig("shap_summary_test_subset.png", dpi=300)
plt.close()

# Full test set SHAP
shap_values_full = explainer.shap_values(X_test_scaled)
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values_full, X_test_scaled, show=False)
plt.title("SHAP Summary Plot (Test Set)")
plt.tight_layout()
plt.savefig("shap_summary_test.png", dpi=300)
plt.close()

# Compute mean SHAP values
if isinstance(shap_values_full, list):
    shap_vals = shap_values_full[1]
else:
    shap_vals = shap_values_full

mean_shap = np.abs(shap_vals).mean(axis=0)
shap_df = pd.DataFrame({
    "CpG": X_test_scaled.columns,
    "Mean_SHAP": mean_shap
}).sort_values(by="Mean_SHAP", ascending=False)

shap_df.to_csv("CpG_SHAP_importance.csv", index=False)
print("SHAP importance saved successfully!")

# -----------------------------
# Full dataset predictions
# -----------------------------
X_scaled_full = pd.DataFrame(scaler.transform(X), columns=X.columns)
y_pred_full = rf.predict(X_scaled_full)

results_rf_full = pd.DataFrame({
    "SampleID": GSE55763_beta["X"],
    "Observed_Age": GSE55763_beta["age"],
    "Pred_RF": y_pred_full
})

results_rf_full.to_csv("RF_Predicted_Ages_full.csv", index=False)
print("Full dataset predictions saved successfully!")
print(results_rf_full.shape)
