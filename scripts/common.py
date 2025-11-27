import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    f1_score,
    accuracy_score,
    precision_score,
    recall_score,
    classification_report
)
from scipy.special import expit # Sigmoid function
import logging

# description: Utility functions for reading fasta sequences, motif extraction, and calculate performance metrics for a classifier.


def fasta_to_dataframe(seq):
    fasta_df = pd.read_csv(seq, sep='>', lineterminator='>', header=None)
    df = fasta_df[0].str.split('\n', n=1, expand=True)
    df[1] = df[1].replace('\n','', regex=True)
    df = df[df[1] != '']
    df = df.dropna()
    df.columns = ['tid','seq']
    return df


def substitute_splice_motifs(seq, positions, dinucleotide, sub='AA'):
    seq = seq.upper()
    seq_list = list(seq)
    for p in positions:
        # Replace motif with AA
        if p + 1 < len(seq_list):
            if (dinucleotide == "GT") & (seq_list[p] == dinucleotide[0]) & (seq_list[p+1] == dinucleotide[1]):
                seq_list[p] = sub[0]
                seq_list[p+1] = sub[1]
            else:
                logging.warning("Position", p, "doesn't match with", dinucleotide)
                
            if (dinucleotide == "AG") & (seq_list[p-1] == dinucleotide[0]) & (seq_list[p] == dinucleotide[1]):
                seq_list[p-1] = sub[0]
                seq_list[p] = sub[1]
            else:
                logging.warning("Position", p, "doesn't match with", dinucleotide)
                
    return ''.join(seq_list)


def get_non_substring_segments(long_str, sub_str):
    start = long_str.find(sub_str)
    if start == -1:
        # substring not found
        return [long_str]
    end = start + len(sub_str)
    # Segments before and after the substring
    segments = []
    if start > 0:
        segments.append(long_str[:start])
    if end < len(long_str):
        segments.append(long_str[end:])
    return segments


def extract_quadruplet_rows(
    df: pd.DataFrame,
    col_name: str,
    first_value: str,
    second_value: str,
    pre_value_any: bool = True,
    pre_specific_value: str = None,
    post_value_any: bool = True,
    post_specific_value: str = None
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    
    """
    Extracts rows forming the pattern X -> first_value -> second_value -> Y
    and returns four separate dataframes: for X, first_value, second_value, and Y.

    Args:
        df (pd.DataFrame): The input DataFrame.
        col_name (str): The name of the column to check.
        first_value (str): The value for the second position in the pattern (e.g., 'A').
        second_value (str): The value for the third position in the pattern (e.g., 'G').
        pre_value_any (bool): If True, X can be any value.
        pre_specific_value (str, optional): The specific value for X if `pre_value_any` is False.
        post_value_any (bool): If True, Y can be any value.
        post_specific_value (str, optional): The specific value for Y if `post_value_any` is False.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            A tuple containing four DataFrames:
            - df_pre: DataFrame for the preceding (X) row.
            - df_first: DataFrame for the 'first_value' (A) row.
            - df_second: DataFrame for the 'second_value' (G) row.
            - df_third: DataFrame for the following (Y) row.
            Returns four empty DataFrames if no such pattern is found,
            or if the input DataFrame is empty, or the column does not exist.
    """
    
    empty_df = pd.DataFrame(columns=df.columns)

    if df.empty:
        print("Input DataFrame is empty.")
        return empty_df, empty_df, empty_df, empty_df

    if col_name not in df.columns:
        raise ValueError(f"Column '{col_name}' not found in the DataFrame.")

    # 1. Identify rows where 'first_value' (A) is present
    is_first_val = (df[col_name] == first_value)

    # 2. Check for 'second_value' (G) immediately after 'first_value' (A)
    is_second_val_after_first = (df[col_name].shift(-1) == second_value)

    # 3. Check for 'third_value' (Y) immediately after 'second_value' (G)
    #    i.e. check the value at current_index + 2
    is_third_val_after_second = pd.Series(False, index=df.index)
    if post_value_any:
        # Check if the value at shift(-2) is not NaN (i.e., exists)
        is_third_val_after_second = ~df[col_name].shift(-2).isna()
    else:
        if post_specific_value is None:
            raise ValueError("`post_specific_value` must be provided if `post_value_any` is False.")
        is_third_val_after_second = (df[col_name].shift(-2) == post_specific_value)

    # Combine these conditions to find the starting index (the 'A' row) of the A->G->Y pattern
    # This identifies the 'first_value' rows that are followed by 'second_value' AND 'third_value'
    first_indices_in_pattern = df.index[is_first_val & is_second_val_after_first & is_third_val_after_second]

    # Filter for the 'pre_value' (X)
    pre_indices_candidate = first_indices_in_pattern - 1
    
    # Filter out pre_indices that are less than 0 (i.e., 'first_value' was at index 0)
    valid_pre_indices = pre_indices_candidate[pre_indices_candidate >= 0]

    # Apply specific value check for 'pre_value' (X) if required
    if not pre_value_any:
        if pre_specific_value is None:
            raise ValueError("`pre_specific_value` must be provided if `pre_value_any` is False.")
        
        actual_pre_values = df.loc[valid_pre_indices, col_name]
        matching_pre_filter = (actual_pre_values == pre_specific_value).values
        valid_pre_indices = valid_pre_indices[matching_pre_filter]

    # Generate the final sets of indices for all four positions
    # These indices correspond to the *start* of each of the 4 elements in the X-A-G-Y pattern.
    final_pre_indices = valid_pre_indices
    final_first_indices = valid_pre_indices + 1
    final_second_indices = valid_pre_indices + 2
    final_third_indices = valid_pre_indices + 3

    # Extract the rows into separate dataframes
    df_pre = df.loc[final_pre_indices].copy()
    df_first = df.loc[final_first_indices].copy()
    df_second = df.loc[final_second_indices].copy()
    df_third = df.loc[final_third_indices].copy()

    return df_pre, df_first, df_second, df_third


def bootstrap_conservation(
    data, 
    genotype, 
    splice_site,
    genotype_col="genotype", 
    splice_site_col="Splice site",
    id_col="id",    
    position_col="position", 
    group_col="group", 
    groups=["Higher usage", "Lower usage"], 
    n_bootstrap=10000
):
    data = data[(data[genotype_col]==genotype) & (data[splice_site_col]==splice_site)].copy()
    
    results = {g: [] for g in groups}
    for i in range(n_bootstrap):
        sample_id = data.drop_duplicates(id_col).sample(frac=0.1, replace=True, random_state=i)[[id_col]]
        sample = pd.merge(sample_id, data)
        for g in groups:
            counts = sample[sample[group_col] == g].value_counts(position_col)
            prop = (counts/sample_id.shape[0]).mean() #.reset_index()
            # prop.columns = [position_col, "proportion"]
            results[g].append(prop)
            
    return results


def performance_metrics(y_true, y_pred, threshold, precision=2):
    """
    Calculates and prints various performance metrics for binary classification.

    Args:
        y_true (array-like): True binary labels.
        y_pred (array-like): Predicted probabilities or scores for the positive class.
        threshold (float): The threshold to convert probabilities to binary labels.
        precision (int, optional): The number of decimal places to format
                                   the float output. Defaults to 4.
    """
    # Create the format specifier string
    format_spec = f".{precision}f"

    # Calculate AUROC and AUPRC
    roc_auc = roc_auc_score(y_true, y_pred)
    auprc = average_precision_score(y_true, y_pred)

    print(f"\nAUROC: {roc_auc:{format_spec}}")
    print(f"AUPRC: {auprc:{format_spec}}")

    # Calculate F1-Score, Accuracy, Precision, Recall
    y_pred_binary_labels = (y_pred >= threshold).astype(int)

    print(f"\nPredicted binary labels unique values: {np.unique(y_pred_binary_labels)}")

    # Calculate metrics
    accuracy = accuracy_score(y_true, y_pred_binary_labels)
    f1 = f1_score(y_true, y_pred_binary_labels, pos_label=1) # Specify pos_label for binary
    precision_val = precision_score(y_true, y_pred_binary_labels, pos_label=1)
    recall = recall_score(y_true, y_pred_binary_labels, pos_label=1)

    print(f"\nMetrics at threshold={threshold}:")
    print(f"Accuracy: {accuracy:{format_spec}}")
    print(f"F1-Score: {f1:{format_spec}}")
    print(f"Precision: {precision_val:{format_spec}}")
    print(f"Recall: {recall:{format_spec}}")

    # Classification Report for a full summary
    print(f"\nClassification Report at threshold={threshold}:")
    print(classification_report(y_true, y_pred_binary_labels, target_names=['Not Positive', 'Positive']))
