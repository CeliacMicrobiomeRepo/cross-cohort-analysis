"""
lodo_xset_train.py

Description:
 - Trials 5 models (SVM, LR, XGBOOST, RF, MLP) with different hyperparameter combinations
 - Performs LODO (Leave-One-Group-Out where groups are datasets) cross-validation with 5 replicates
 - Determines the best model and hyperparameters based on mean AUC across folds and replicates
 - Optionally performs cross-dataset validation (XSET) using the best hyperparameters of each of the 5 models

Parts:
 1. Hyperparameter grid search 
 2. Write model_results.txt (best mean AUC and hyperparameter combination for each 5 models)
 3. Best model testing (same as 1, but on the best hyperparameters only)
 4. Write best_model_results.txt (best mean AUC for each of the 5 model's best hyperparameters)
 5. Write ROC curve plots (.png) for each of the 5 models
 6. Write ROC curve data (.tsv) for each of the 5 models
 7. Write summary_best_model.tsv (a tiny tsv file with 3 values: auc, model_name, hyperparameters)
 8. Write about_train.txt (a text file with general info about the training e.g. date, number of replicates, etc.)
 9. Write Comprehensive results for every hyperparameter combination for every model (comprehensive_hyperparameter_results.tsv)
 10. Cross-dataset validation (OPTIONAL)
    a) Perform cross-dataset validation (XSET) (For all 5 models using the best hyperparameters...
        trains and evaluates on all pairs of datasets including self-testing). This results in...
        a matrix of mean AUC values for each of the 5 model's best hyperparameters.
    b) Write xset_results_{MODEL}.tsv for each of the 5 models. (rows are training datasets, columns are testing datasets, values are mean AUCs)
    c) Write mean_xset_results.tsv (The mean of each of the 5 model's mean AUCs across all cross-dataset combinations (excluding self-validation))

For input, this script takes the output of either:
    - construct_dataset_ASV.py 
    - construct_dataset_taxa.py

This script is designed to be run from the command line inside the dataset directory.

Example:
```
conda activate all_env
cd /home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_main/duodenum_active_log10_after/
python ../../03_lodo_xset_train.py
```

"""


# Imports ------------------------------------------
import pandas as pd
import numpy as np
import datetime
import os
from sklearn.model_selection import LeaveOneGroupOut, ParameterGrid
from sklearn.metrics import roc_auc_score
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from custom_sklearn import MLPDropout
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
import time

# Record the start time
start_time = time.time()


# Options ------------------------------------------
# Path to the dataset directory (also output directory)
DATASET_DIR_PATH = os.getcwd()   # e.g. "./machine_learning/datasets_main/duodenum_active_log10_after/"
# Extract the analysis group from the dataset directory
# e.g.   'duodenum_active'   'stool_prospective'   'stool_active'   'stool_treated'
GROUP_NAME = "_".join(DATASET_DIR_PATH.split("/")[-1].split("_")[0:2])
# Path to the feature matrix TSV file (input)
FEATURES_TSV_PATH = os.path.join(DATASET_DIR_PATH, "sample_asv_abundances.tsv")
# Path to the labels TSV file (input)
LABELS_TSV_PATH = os.path.join(DATASET_DIR_PATH, "sample_labels.tsv")
# Output subdirectory
OUTPUT_SUBDIR_PATH = os.path.join(DATASET_DIR_PATH, "lodo_results/")
# Create the output subdirectory if it doesn't exist
if not os.path.exists(OUTPUT_SUBDIR_PATH):
    os.makedirs(OUTPUT_SUBDIR_PATH)

# Number of replicates for cross-validation
NUM_REPLICATES = 5

# Hyperparameter grids for each model
hyperparameter_grids = {
    'mlp': {
        'hidden_layer_sizes': [(64,), (256,), (128, 64)],
        'dropout': [None, 0.25],
        'learning_rate_init': [0.01, 0.001],
        'learning_rate': ['adaptive'],
        'max_iter': [2000],
        'batch_size': ['auto']
    },
    'xgboost': {
        'n_estimators': [100, 500],
        'learning_rate': [0.1, 0.01],
        'gamma': [0, 0.5],
        'colsample_bytree': [0.8, 1]
    },
    'svm': {
        'C': [0.1, 1],
        'kernel': ['linear', 'rbf', 'sigmoid'],
        'gamma': ['scale', 0.1, 1]
    },
    'lr': {
        'C': [0.01, 0.1, 1],
        'penalty': ['l1', 'l2'],
        'solver': ['lbfgs', 'saga', 'liblinear', 'newton-cg'],
        'max_iter': [10000]
    },
    'rf': {
        'n_estimators': [200, 750],
        'max_depth': [15, 25],
        'max_features': ["sqrt", "log2"]
    }
}

# Mapping from model names to model classes
model_classes = {
    'mlp': MLPDropout,
    'xgboost': XGBClassifier,
    'svm': SVC,
    'lr': LogisticRegression,
    'rf': RandomForestClassifier,
}

# Do cross dataset validation?
DOING_CROSS_DATASET_VALIDATION = True



# Setup ------------------------------------------

# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_SUBDIR_PATH):
    os.makedirs(OUTPUT_SUBDIR_PATH)

# Load feature matrix and labels
X = pd.read_csv(FEATURES_TSV_PATH, sep='\t', index_col=0)
y = pd.read_csv(LABELS_TSV_PATH, sep='\t')

# Ensure that the sample IDs match between features and labels
X = X.sort_index()
y = y.set_index('Sample_ID')
y = y.loc[X.index]

# Extract labels and group IDs for cross-validation
labels = y['Disease_Status'].astype(str).values  # Convert booleans to strings
groups = y['Dataset_ID'].values  # For leave-one-dataset-out cross-validation

# Encode labels as integers
le = LabelEncoder()
y_encoded = le.fit_transform(labels)

# Clean up column names to remove unsupported characters for XGBoost
X.columns = [str(col).replace('[', '').replace(']', '').replace('<', '').replace('>', '') for col in X.columns]

# Get unique datasets
unique_datasets = np.unique(groups)





# Hyperparameter grid search ------------------------------------------
# For every hyperparameter combination of each model, perform LODO cross-validation NUM_REPLICATES times
# Store the mean AUC for each hyperparameter combination -> this tells us the best hyperparameters for each model

# To store results for each model
model_results = {}

# To store best model information
best_model_name = None
best_model_mean_auc = -np.inf
best_model_results = None

# Perform model selection
for model_name in ['mlp', 'xgboost', 'svm', 'lr', 'rf']:
    print(f"\nEvaluating model: {model_name}")
    ModelClass = model_classes[model_name]
    param_grid = hyperparameter_grids[model_name]
    param_combinations = list(ParameterGrid(param_grid))
    best_hyperparams = None
    best_mean_auc = -np.inf
    best_replicate_aucs = []
    # For storing mean AUC across replicates for each hyperparameter combination
    hyperparam_auc_results = []
    for hyperparams in param_combinations:
        # Create distinct seeds for each replicate
        replicate_rng = np.random.RandomState()
        replicate_seeds = [replicate_rng.randint(0, 1e6) for _ in range(NUM_REPLICATES)]
        replicate_mean_aucs = []
        for replicate in range(1, NUM_REPLICATES +1):
            # Set random seed for this replicate
            base_seed = replicate_seeds[replicate - 1]
            MODEL_SEED = base_seed + 1
            DATA_SEED = base_seed + 2
            WEIGHT_INIT_SEED = base_seed + 3
            # Shuffle data
            rng = np.random.RandomState(DATA_SEED)
            shuffled_indices = rng.permutation(len(X))
            X_shuffled = X.iloc[shuffled_indices].reset_index(drop=True)
            y_shuffled = y_encoded[shuffled_indices]
            groups_shuffled = groups[shuffled_indices]
            # Initialize the model with the given hyperparameters
            if model_name in ['rf', 'xgboost']:
                model = ModelClass(**hyperparams, n_jobs=1)
            else:
                model = ModelClass(**hyperparams)
            if model_name == 'mlp':
                # Add specific MLP initialization randomization
                model.set_params(random_state=MODEL_SEED)
                model._random_state = WEIGHT_INIT_SEED
            if hasattr(model, 'random_state'):
                model.set_params(random_state=MODEL_SEED)
            if hasattr(model, 'seed'):
                model.set_params(seed=MODEL_SEED)
            if hasattr(model, 'probability'):
                model.set_params(probability=True)
            # Initialize leave-one-dataset-out cross-validator
            lodo = LeaveOneGroupOut()
            # Store AUCs for this replicate
            replicate_auc_scores = []
            # Perform cross-validation
            for train_idx, test_idx in lodo.split(X_shuffled, y_shuffled, groups_shuffled):
                X_train, X_test = X_shuffled.iloc[train_idx], X_shuffled.iloc[test_idx]
                y_train, y_test = y_shuffled[train_idx], y_shuffled[test_idx]
                try:
                    # Fit the model
                    model.fit(X_train, y_train)
                except ValueError as e:
                    replicate_auc_scores.append(0)  # Set AUC to 0 for this hyperparameter combination
                    continue
                # Make predictions
                if hasattr(model, "predict_proba"):
                    y_pred_proba = model.predict_proba(X_test)[:, 1]
                elif hasattr(model, "decision_function"):
                    y_pred_proba = model.decision_function(X_test)
                else:
                    # Cannot compute AUC without probabilities or decision function
                    continue
                # Evaluate performance
                auc = roc_auc_score(y_test, y_pred_proba)
                # Append AUC
                replicate_auc_scores.append(auc)
            # Compute mean AUC across groups for this replicate
            if len(replicate_auc_scores) > 0:
                mean_auc = np.mean(replicate_auc_scores)
                replicate_mean_aucs.append(mean_auc)
            else:
                mean_auc = np.nan
                replicate_mean_aucs.append(mean_auc)
        # Compute mean AUC across replicates for this hyperparameter combination
        if len(replicate_mean_aucs) > 0:
            hyperparam_mean_auc = np.nanmean(replicate_mean_aucs)
        else:
            hyperparam_mean_auc = np.nan
        hyperparam_auc_results.append({
            'hyperparameters': hyperparams,
            'mean_auc': hyperparam_mean_auc
        })
        # Update best hyperparameters for this model
        if hyperparam_mean_auc > best_mean_auc:
            best_mean_auc = hyperparam_mean_auc
            best_hyperparams = hyperparams
            best_replicate_aucs = replicate_mean_aucs.copy()
        print(f"  AUC: {hyperparam_mean_auc:.3f}  For hyperparameters: {hyperparams}")
    # Store best hyperparameters and mean AUC for this model
    model_results[model_name] = {
        'best_hyperparameters': best_hyperparams,
        'mean_auc': best_mean_auc,
        'replicate_aucs': best_replicate_aucs,
        'hyperparam_auc_results': hyperparam_auc_results
    }
    # Update best model overall
    if best_mean_auc > best_model_mean_auc:
        best_model_mean_auc = best_mean_auc
        best_model_name = model_name
        best_model_results = {
            'replicate_aucs': best_replicate_aucs,
            'best_hyperparameters': best_hyperparams
        }




# Write model_results.txt ------------------------------------------
model_results_txt_path = os.path.join(OUTPUT_SUBDIR_PATH, 'model_results.txt')
with open(model_results_txt_path, 'w') as f:
    # Sort model names alphabetically
    for model_name in sorted(model_results.keys()):
        result = model_results[model_name]
        f.write(f"Model: {model_name}\n")
        f.write(f"Best Hyperparameters: {result['best_hyperparameters']}\n")
        f.write(f"Mean AUC across replicates: {result['mean_auc']:.4f}\n\n")

print(f"Model results saved to {model_results_txt_path}")




# Collect Best Model Results ----------------------------------
# For each model, perform LODO cross-validation NUM_REPLICATES times on the best hyperparameters
# We run these again so that we know the high performance was not due to random chance
# Store the mean AUC for each model -> this tells us the performance of the best hyperparameters of the 5 models
# We also collect a bunch of useful data and objects for later analysis

# Prepare to store per-replicate and per-group AUCs for all best models
all_best_models_results = []

# To store results for each best model
best_model_results = {}

# To store ROC curves
roc_curves = {}

# Initialize variables for tracking the best model overall
best_model_name = None
best_hyperparams = None
best_model_mean_auc = -np.inf

for model_name in model_classes.keys():
    print(f"\nProcessing best model results for: {model_name}")
    best_model_class = model_classes[model_name]
    best_hyperparams = model_results[model_name]['best_hyperparameters']
    per_replicate_group_aucs = []

    # Create distinct seeds for each replicate
    replicate_rng = np.random.RandomState()
    replicate_seeds = [replicate_rng.randint(0, 1e6) for _ in range(NUM_REPLICATES)]

    for replicate in range(1, NUM_REPLICATES + 1):
        # Set random seeds for this replicate
        base_seed = replicate_seeds[replicate - 1]
        MODEL_SEED = base_seed + 1
        DATA_SEED = base_seed + 2
        WEIGHT_INIT_SEED = base_seed + 3
        # Shuffle data
        rng = np.random.RandomState(DATA_SEED)
        shuffled_indices = rng.permutation(len(X))
        X_shuffled = X.iloc[shuffled_indices].reset_index(drop=True)
        y_shuffled = y_encoded[shuffled_indices]
        groups_shuffled = groups[shuffled_indices]
        # Initialize the model with the best hyperparameters
        model = best_model_class(**best_hyperparams)
        if model_name == 'mlp':
            # Add specific MLP initialization randomization
            model.set_params(random_state=MODEL_SEED)
            model._random_state = WEIGHT_INIT_SEED
        if hasattr(model, 'random_state'):
            model.set_params(random_state=MODEL_SEED)
        if hasattr(model, 'seed'):
            model.set_params(seed=MODEL_SEED)
        if hasattr(model, 'probability'):
            model.set_params(probability=True)
        # Initialize leave-one-dataset-out cross-validator
        lodo = LeaveOneGroupOut()
        # Store AUCs for this replicate per group
        replicate_group_aucs = {}
        # Perform cross-validation
        for train_idx, test_idx in lodo.split(X_shuffled, y_shuffled, groups_shuffled):
            X_train, X_test = X_shuffled.iloc[train_idx], X_shuffled.iloc[test_idx]
            y_train, y_test = y_shuffled[train_idx], y_shuffled[test_idx]
            test_group = groups_shuffled[test_idx][0]
            # Fit the model
            model.fit(X_train, y_train)
            # Make predictions
            if hasattr(model, "predict_proba"):
                y_pred_proba = model.predict_proba(X_test)[:, 1]
            elif hasattr(model, "decision_function"):
                y_pred_proba = model.decision_function(X_test)
            else:
                # Cannot compute AUC without probabilities or decision function
                continue
            # Store ROC curve for this fold
            if model_name not in roc_curves:
                roc_curves[model_name] = {'fpr': [], 'tpr': []}
            fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
            roc_curves[model_name]['fpr'].append(fpr)
            roc_curves[model_name]['tpr'].append(tpr)
            # Evaluate performance
            auc = roc_auc_score(y_test, y_pred_proba)
            # Store AUC for this group
            replicate_group_aucs[test_group] = auc
        # Compute mean AUC across groups for this replicate
        replicate_mean_auc = np.mean(list(replicate_group_aucs.values()))
        replicate_group_aucs['Mean'] = replicate_mean_auc
        replicate_group_aucs['Model'] = model_name
        per_replicate_group_aucs.append(replicate_group_aucs)

    # Create a DataFrame to store the results for this model
    model_results_df = pd.DataFrame(per_replicate_group_aucs)
    # Add a column for replicate number, explicitly as string type
    model_results_df['Replicate'] = range(1, len(model_results_df) + 1)
    model_results_df['Replicate'] = model_results_df['Replicate'].astype(str)

    # Compute average over replicates
    average_row = model_results_df.mean(numeric_only=True)
    average_row['Replicate'] = 'Average'
    average_row['Model'] = model_name
    average_row = pd.DataFrame([average_row])
    model_results_df = pd.concat([model_results_df, average_row], ignore_index=True)

    # Append results to the overall list
    all_best_models_results.append(model_results_df)

    # Check if this model's mean AUC is the best overall
    model_mean_auc = average_row['Mean'].item()
    if model_mean_auc > best_model_mean_auc:
        best_model_mean_auc = model_mean_auc
        best_model_name = model_name
        best_hyperparams = best_hyperparams

    # Store best model results
    best_model_results[model_name] = {
        'auc': model_mean_auc,
        'hyperparameters': best_hyperparams
    }

# Combine results for all models into a single DataFrame
all_best_models_df = pd.concat(all_best_models_results, ignore_index=True)

# Reorder columns to put 'Replicate' first
all_best_models_df = all_best_models_df[['Replicate'] + [col for col in all_best_models_df.columns if col != 'Replicate']]
# Save to TSV
all_best_models_results_tsv_path = os.path.join(OUTPUT_SUBDIR_PATH, 'best_models_results.tsv')
all_best_models_df.to_csv(all_best_models_results_tsv_path, sep='\t', index=False)

print(f"All best model results saved to {all_best_models_results_tsv_path}")




# Write best_model_results.txt ----------------------------------
best_model_results_txt_path = os.path.join(OUTPUT_SUBDIR_PATH, 'best_model_results.txt')
with open(best_model_results_txt_path, 'w') as f:
    # Sort model names alphabetically
    for model_name in sorted(best_model_results.keys()):
        result = best_model_results[model_name]
        f.write(f"Model: {model_name}\n")
        f.write(f"Best Hyperparameters: {result['hyperparameters']}\n")
        f.write(f"Mean AUC across replicates: {result['auc']:.4f}\n\n")

print(f"Model results saved to {model_results_txt_path}")



# Write best models ROC curves --------------------------------
# Create the best models subdirectory if it doesn't exist
BEST_MODELS_SUBDIR_PATH = os.path.join(OUTPUT_SUBDIR_PATH, "best_models/")
if not os.path.exists(BEST_MODELS_SUBDIR_PATH):
    os.makedirs(BEST_MODELS_SUBDIR_PATH)

# Write ROC curves for best models
plt.style.use('seaborn-v0_8-colorblind')
for model_name in roc_curves:
    plt.figure(figsize=(8, 6))
    
    # Plot individual ROC curves (light color)
    for fpr, tpr in zip(roc_curves[model_name]['fpr'], roc_curves[model_name]['tpr']):
        plt.plot(fpr, tpr, alpha=0.1, color='lightblue')
    
    # Calculate and plot mean ROC curve (bold)
    # Interpolate all ROC curves to same x-axis points
    mean_fpr = np.linspace(0, 1, 100)
    tprs = []
    for fpr, tpr in zip(roc_curves[model_name]['fpr'], roc_curves[model_name]['tpr']):
        tprs.append(np.interp(mean_fpr, fpr, tpr))
    mean_tpr = np.mean(tprs, axis=0)
    
    plt.plot(mean_fpr, mean_tpr, 'b-', 
             label=f'{model_name.upper()} (Mean AUC = {best_model_results[model_name]["auc"]:.3f})',
             linewidth=2)
    
    # Add diagonal line
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    
    # Customize plot
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curves for {model_name.upper()}')
    plt.legend(loc='lower right')
    
    # Save plot
    roc_path = os.path.join(BEST_MODELS_SUBDIR_PATH, f'roc_{model_name}.png')
    plt.savefig(roc_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"ROC curve for {model_name} saved to {roc_path}")



# Write best models ROC data --------------------------------
for model_name in roc_curves:
    # We'll use the same mean_fpr and mean_tpr we calculated for the bold (mean) ROC curve
    mean_fpr = np.linspace(0, 1, 100)
    tprs = []
    for fpr_vals, tpr_vals in zip(roc_curves[model_name]['fpr'], roc_curves[model_name]['tpr']):
        tprs.append(np.interp(mean_fpr, fpr_vals, tpr_vals))
    mean_tpr = np.mean(tprs, axis=0)

    # Save the mean ROC curve data to TSV
    roc_data_df = pd.DataFrame({
        'fpr': mean_fpr,
        'tpr': mean_tpr
    })
    roc_data_path = os.path.join(BEST_MODELS_SUBDIR_PATH, f'roc_{model_name}_data.tsv')
    roc_data_df.to_csv(roc_data_path, sep='\t', index=False)
    print(f"ROC curve data for {model_name} saved to {roc_data_path}")



# Write summary_best_model.tsv ----------------------------------
summary_best_model_tsv_path = os.path.join(OUTPUT_SUBDIR_PATH, 'summary_best_model.tsv')
# Simply three cells: auc, model_name, hyperparameters
text = f"{best_model_mean_auc:.4f}\t{best_model_name.upper()}\t{best_model_results[best_model_name]['hyperparameters']}"
with open(summary_best_model_tsv_path, 'w') as f:
    f.write(text)
print(f"Summary best model saved to {summary_best_model_tsv_path}")



# Write about_train.txt ----------------------------------
with open(os.path.join(OUTPUT_SUBDIR_PATH, "about_train.txt"), "w") as about_file:
    # Date
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    about_file.write(f"Training Summary\n")
    about_file.write(f"Date and Time: {current_date}\n\n")
    
    # Input File Paths
    about_file.write("Input Files:\n")
    about_file.write(f"  - Feature Matrix: {FEATURES_TSV_PATH}\n")
    about_file.write(f"  - Labels File: {LABELS_TSV_PATH}\n\n")
    
    # Cross-Validation and Replicates
    about_file.write("Training Configuration:\n")
    about_file.write(f"  - Cross-validation Method: Leave-One-Dataset-Out\n")
    about_file.write(f"  - Number of Replicates: {NUM_REPLICATES}\n\n")
    
    # Hyperparameter Grids
    about_file.write("Hyperparameter Grids:\n")
    for model_name, params in hyperparameter_grids.items():
        about_file.write(f"  {model_name.upper()}:\n")
        for param, values in params.items():
            about_file.write(f"    - {param}: {values}\n")
        about_file.write("\n")
    
    # Models Used
    about_file.write("Models Evaluated:\n")
    for model_name in model_classes.keys():
        about_file.write(f"  - {model_name.upper()}\n")
    
    # Additional Notes
    about_file.write("\nNotes:\n")
    about_file.write("  - This training session uses randomized replicates with reshuffled data.\n")
    about_file.write("  - Best model and hyperparameters are selected based on mean AUC across replicates.\n")



# Write comprehensive_hyperparameter_results.tsv ----------------------------------
all_results = []
for model_name, hyperparam_auc_results in model_results.items():
    for result in hyperparam_auc_results['hyperparam_auc_results']:
        # Append each result as a row
        row = {
            'model': model_name,
            **result['hyperparameters'],  # Unpack hyperparameters into columns
            'mean_auc': result['mean_auc']
        }
        all_results.append(row)

# Create a DataFrame from the results
results_df = pd.DataFrame(all_results)

# Reorder columns to ensure 'model' and 'mean_auc' are first
columns_order = ['model', 'mean_auc'] + [col for col in results_df.columns if col not in ['model', 'mean_auc']]
results_df = results_df[columns_order]

# Save the results to a TSV file
results_tsv_path = os.path.join(OUTPUT_SUBDIR_PATH, 'comprehensive_hyperparameter_results.tsv')
results_df.to_csv(results_tsv_path, sep='\t', index=False)

print(f"All model hyperparameter results saved to {results_tsv_path}")

# Total elapsed time
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total time taken: {elapsed_time:.2f} seconds")



# Summary ------------------
print(f"\nBest Model: {best_model_name}")
print(f"Best Hyperparameters: {best_hyperparams}")
print(f"Mean AUC across replicates: {best_model_mean_auc:.4f}")




# Cross dataset validation ----------------------------------
# Pre-generate replicate seeds for cross-dataset validation
replicate_seeds = np.random.randint(0, 1e6, size=NUM_REPLICATES)
if DOING_CROSS_DATASET_VALIDATION:
    print("\nStarting cross-dataset validation...")

    # This stores the mean AUCs for each model (except when testing on the same dataset as training)
    all_non_self_val_aucs = {
        'mlp': [],
        'xgboost': [],
        'svm': [],
        'lr': [],
        'rf': []
    }

    for model_name in ['mlp', 'xgboost', 'svm', 'lr', 'rf']:
        best_hyperparams = model_results[model_name]['best_hyperparameters']
        ModelClass = model_classes[model_name]
        xset_matrix = np.zeros((len(unique_datasets), len(unique_datasets)))

        for i, train_dataset in enumerate(unique_datasets):
            # Create train set
            train_mask = (groups == train_dataset)
            X_train = X[train_mask].reset_index(drop=True)
            y_train = y_encoded[train_mask]

            for j, test_dataset in enumerate(unique_datasets):
                test_mask = (groups == test_dataset)
                X_test = X[test_mask].reset_index(drop=True)
                y_test = y_encoded[test_mask]

                # Repeat to get mean AUC
                replicate_aucs = []
                for replicate_idx in range(NUM_REPLICATES):
                    # Set random seeds for this replicate
                    base_seed = replicate_seeds[replicate_idx]
                    MODEL_SEED = base_seed + 1
                    DATA_SEED = base_seed + 2
                    WEIGHT_INIT_SEED = base_seed + 3
                    rng = np.random.RandomState(DATA_SEED)
                    shuffled_idx = rng.permutation(len(X_train))
                    X_train_shuffled = X_train.iloc[shuffled_idx].reset_index(drop=True)
                    y_train_shuffled = y_train[shuffled_idx]

                    # Build model and fit
                    if model_name in ['rf', 'xgboost']:
                        model = ModelClass(**best_hyperparams, n_jobs=1)
                    else:
                        model = ModelClass(**best_hyperparams)
                    if model_name == 'mlp':
                        # MLP seeds
                        model.set_params(random_state=MODEL_SEED)
                        model._random_state = WEIGHT_INIT_SEED
                    if hasattr(model, 'random_state'):
                        model.set_params(random_state=MODEL_SEED)
                    if hasattr(model, 'seed'):
                        model.set_params(seed=MODEL_SEED)
                    if hasattr(model, 'probability'):
                        model.set_params(probability=True)

                    model.fit(X_train_shuffled, y_train_shuffled)

                    # Predict probabilities and compute AUC
                    if hasattr(model, 'predict_proba'):
                        y_pred_proba = model.predict_proba(X_test)[:, 1]
                    elif hasattr(model, 'decision_function'):
                        y_pred_proba = model.decision_function(X_test)
                    else:
                        # Can't compute AUC otherwise
                        continue

                    replicate_aucs.append(roc_auc_score(y_test, y_pred_proba))
                # Get the mean AUC for all replicates
                mean_of_replicates = np.mean(replicate_aucs)
                # Store the mean AUC for all replicates
                xset_matrix[i, j] = mean_of_replicates
                # Unless the test and train datasets are the same
                if i != j:
                    # Store the mean AUC for all replicates (to get mean later on)
                    all_non_self_val_aucs[model_name].append(mean_of_replicates)

        # Convert xset_matrix to a DataFrame, write to file
        xset_df = pd.DataFrame(xset_matrix, index=unique_datasets, columns=unique_datasets)
        xset_path = os.path.join(OUTPUT_SUBDIR_PATH, f"xset_results_{model_name}.tsv")
        xset_df.to_csv(xset_path, sep='\t')
        print(f"Cross-dataset results for {model_name} saved to {xset_path}")

    # Compute the mean AUC for each model (Excluding self-validation) and save to TSV
    mean_aucs = []
    for model_name in all_non_self_val_aucs:
        mean_auc = np.mean(all_non_self_val_aucs[model_name])
        mean_aucs.append({
            'Model': model_name,
            'Mean AUC': mean_auc
        })
        print(f"Mean AUC for {model_name} (excluding self-validation): {mean_auc:.4f}")
    # Create dataframe and sort by Mean AUC
    mean_aucs_df = pd.DataFrame(mean_aucs)
    mean_aucs_df = mean_aucs_df.sort_values('Mean AUC', ascending=False)
    # Save to TSV file
    mean_xset_path = os.path.join(OUTPUT_SUBDIR_PATH, f"mean_xset_results.tsv")
    mean_aucs_df.to_csv(mean_xset_path, sep='\t', index=False)
    print(f"Mean cross-dataset results saved to {mean_xset_path}")


# Total elapsed time
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total time taken: {elapsed_time:.2f} seconds")

