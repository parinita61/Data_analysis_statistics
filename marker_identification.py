# Import necessary libraries
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    precision_score, 
    recall_score, 
    roc_auc_score, 
    confusion_matrix, 
    roc_curve, 
    auc
)
from sklearn.preprocessing import LabelEncoder
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import Parallel, delayed  # Import joblib for parallel execution


# Load the dataset
data = pd.read_csv("PupilBioTest_PMP_revA.csv")

# Shuffle the data to ensure random distribution of tissue types
data = data.sample(frac=1, random_state=42).reset_index(drop=True)

# Encode categorical columns
label_encoder = LabelEncoder()
data['strand_encoded'] = label_encoder.fit_transform(data['strand'])
data['Tissue_encoded'] = label_encoder.fit_transform(data['Tissue'])

# Extract features (CpG Coordinates and methylation patterns) and labels (Tissue type)
features = ['strand_encoded', 'CpG_Coordinates', '`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
X = data[features]
y = data['Tissue_encoded']

# Split data into training and testing sets (randomly shuffled)
X_train, X_test, y_train, y_test = train_test_split(
    X.drop(columns=['CpG_Coordinates']),  # Drop non-predictive column
    y, 
    test_size=0.3, 
    stratify=y, 
    random_state=42
)

# Train a Random Forest Classifier
model = RandomForestClassifier(random_state=42, n_jobs=-1)  # n_jobs=-1 uses all available cores
model.fit(X_train, y_train)

# Predict labels and probabilities
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)

# Calculate evaluation metrics
precision = precision_score(y_test, y_pred, average='macro')
recall = recall_score(y_test, y_pred, average='macro')

try:
    # Check if y_test is a 1D array, and flatten y_prob to 1D if necessary
    if y_prob.shape[1] > 1:  # More than one class
        # Compute AUC for multiclass using the "one-vs-rest" (ovr) approach
        roc_auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')
    else:
        # If only binary classification, directly compute AUC
        roc_auc = roc_auc_score(y_test, y_prob)
except Exception as e:
    print(f"Error computing AUC: {e}")
#roc_auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')


# Assign statistical confidence (e.g., Fisher's exact test) to each PMP

# Define the Fisher's Exact Test function
def fisher_test(row):
    try:
        # Create the contingency table for the Fisher's Exact Test
        contingency_table = np.array([[row['`000'], row['`111']], [row['`001'], row['`110']]])
        if contingency_table.shape == (2, 2):  # Ensure it's a 2x2 table
            return fisher_exact(contingency_table, alternative='two-sided')[1]
        else:
            return np.nan  # Return NaN if the table is not 2x2
    except Exception as e:
        return np.nan  # Return NaN if there's an exception

# Parallelize the Fisher's Exact Test across the rows of the dataset
p_values = Parallel(n_jobs=-1)(delayed(fisher_test)(row) for _, row in data.iterrows())
data['p_value'] = p_values

# Filter PMPs with high specificity for tissue differentiation
significant_pmp = data[data['p_value'] < 0.05]

# Calculate the mean variant read fraction (VRF) for each PMP in both tissues
vrf_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
data['VRF'] = data[vrf_columns].sum(axis=1)
mean_vrf = data.groupby('Tissue')['VRF'].mean().reset_index()
mean_vrf.columns = ['Tissue', 'Mean_VRF']

# Save results
significant_pmp.to_csv("significant_pmps.csv", index=False)
mean_vrf.to_csv("mean_vrf.csv", index=False)

# Display metrics and key results
print(f"Precision: {precision:.2f}")
print(f"Recall: {recall:.2f}")
print(f"ROC AUC: {roc_auc:.2f}")
print("\nSignificant PMPs saved to 'significant_pmps.csv'")
print("Mean VRF for each tissue saved to 'mean_vrf.csv'")

# --- Plotting Section ---

# 1. ROC Curve
plt.figure(figsize=(8, 6))
for i, class_label in enumerate(label_encoder.classes_):
    fpr, tpr, _ = roc_curve((y_test == i).astype(int), y_prob[:, i])  # One-vs-rest
    roc_auc_value = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=2, label=f'Class {class_label} (AUC = {roc_auc_value:.2f})')

plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc='lower right')
plt.savefig("roc_curve.png")  # Save the ROC curve plot
plt.show()

# 2. Confusion Matrix Plot
conf_matrix = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(6, 5))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', cbar=False, 
            xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title('Confusion Matrix')
plt.savefig("confusion_matrix.png")  # Save the confusion matrix plot
plt.show()

# 3. Feature Importance Plot
feature_importances = model.feature_importances_
sorted_idx = np.argsort(feature_importances)[::-1]

plt.figure(figsize=(10, 6))
plt.barh(range(len(sorted_idx)), feature_importances[sorted_idx], align='center')
plt.yticks(range(len(sorted_idx)), np.array(X_train.columns)[sorted_idx])
plt.xlabel('Feature Importance')
plt.title('Feature Importance Plot')
plt.savefig("feature_importance.png")  # Save the feature importance plot
plt.show()
