from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import normalized_mutual_info_score
import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import logging

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')



def nonss(fasta, chrom, start, end):
    fd = fasta[fasta.tid==chrom].copy()
    sequence = fd.seq.iloc[0].upper()
    
    positions = []
    for i in range(int(start), int(end)):
        if sequence[i:i+2] == "GT":
            positions.append([chrom, i, i+2, "GT(non-donor)", "none", "+"])
        elif sequence[i:i+2] == "AG":
            positions.append([chrom, i, i+2, "AG(non-acceptor)", "none", "+"])
        elif sequence[i:i+2] == "AC":
            positions.append([chrom, i, i+2, "GT(non-donor)", "none", "-"])
        elif sequence[i:i+2] == "CT":
            positions.append([chrom, i, i+2, "AG(non-acceptor)", "none", "-"])
            
    return pd.DataFrame(positions)



def cal_metric_by_group(labels, preds, metric_fun=normalized_mutual_info_score, donors: list=None, acceptors: list=None, by_group: bool=True):
    """
    Calculates a metric (e.g., NMI) for specified groups (donors/acceptors) or on all labels.

    Args:
        labels (array-like): True labels for the data points.
        preds (array-like): Predicted labels/clusters for the data points.
        metric_fun (callable, optional): The metric function to apply.
                                        Defaults to normalized_mutual_info_score.
                                        It should accept two array-like arguments (y_true, y_pred).
        donors (list, optional): A list of labels to be considered "donors".
                                If provided, metric will be calculated for these labels.
        acceptors (list, optional): A list of labels to be considered "acceptors".
                                    If provided, metric will be calculated for these labels.
        by_group (bool, optional): If True, metrics are calculated for 'donors' and/or
                                   'acceptors' groups. If False, the metric is calculated
                                   on all labels. Defaults to True.

    Returns:
        tuple or float:
            - If by_group is True and both donors and acceptors provided:
              (average_score, donor_score, acceptor_score)
            - If by_group is True and only one of donors/acceptors provided:
              the score for that group.
            - If by_group is False: the score for all labels.

    Raises:
        ValueError: If by_group is True but no valid donor or acceptor labels are found.
    """
    if by_group:
        score1 = None
        score2 = None

        if donors:
            k1 = np.isin(labels, donors)
            if np.any(k1): # Ensure there are actual donor labels to calculate score
                score1 = metric_fun(labels[k1], preds[k1])
            else:
                logging.warning("'donors' list provided but no matching labels found.")


        if acceptors:
            k2 = np.isin(labels, acceptors)
            if np.any(k2): # Ensure there are actual acceptor labels to calculate score
                score2 = metric_fun(labels[k2], preds[k2])
            else:
                logging.warning("'acceptors' list provided but no matching labels found.")

        if score1 is not None and score2 is not None:
            return (score1 + score2) / 2, score1, score2
        elif score1 is not None:
            return score1
        elif score2 is not None:
            return score2
        else:
            # Both donors and acceptors were provided but no matching labels,
            # or neither list was provided and by_group was True.
            raise ValueError("No valid donor or acceptor labels found for metric calculation.")

    else:
        # Calculate metric on all labels
        return metric_fun(labels, preds)
    


def seq_to_dnabert_kmers(seq, k: int):
    kmers = list()
    for i in range(0, len(seq) - k + 1):
        kmers.append(seq[i:i+k])
    return ' '.join(kmers)



def prepare_single_sequence_input(
    sequence: str,
    tokenizer: AutoTokenizer,
    max_len: int = 400,
    dnabert_k: int = None
) -> dict:
    """
    Prepares a single DNA sequence for input to a SpliceBERT model.

    Args:
        sequence (str): The DNA sequence string (e.g., 'AGCT...').
        tokenizer (AutoTokenizer): The tokenizer loaded from the SpliceBERT model.
        max_len (int): The maximum sequence length the model expects.
        dnabert_k (int, optional): K-mer size for DNA-BERT specific tokenization.
                                   If None, assumes 1-mer (single base) tokenization.

    Returns:
        dict: A dictionary containing 'input_ids' and 'attention_mask' tensors,
              ready for model inference.
    """
    # 1. Uppercase the sequence
    seq_upper = sequence.upper()

    # 2. Handle sequence length and cropping (mimicking SpliceatorDataset)
    if len(seq_upper) > max_len:
        # Crop from the middle to max_len
        skip_left = (len(seq_upper) - max_len) // 2
        processed_seq = seq_upper[skip_left:skip_left + max_len]
    else:
        # If sequence is shorter than max_len, the tokenizer will handle padding
        processed_seq = seq_upper

    # 3. K-merize if dnabert_k is provided
    if dnabert_k is not None:
        processed_seq_for_tokenizer = seq_to_dnabert_kmers(processed_seq, k=dnabert_k)
    else:
        # For 1-mer tokenization, tokenizer expects space-separated characters
        processed_seq_for_tokenizer = ' '.join(list(processed_seq))

    # 4. Tokenization
    # `return_tensors='pt'` makes it return PyTorch tensors
    # `padding='max_length'` pads to max_len if shorter
    # `truncation=True` truncates if longer (though we pre-cropped)
    encoding = tokenizer(
        processed_seq_for_tokenizer,
        return_tensors='pt',
        padding='max_length', # Pad to max_len if the sequence is shorter
        truncation=True,      # Truncate if it's still longer (shouldn't happen with pre-cropping)
        max_length=max_len    # Ensure consistent length
    )

    return {
        'input_ids': encoding['input_ids'].squeeze(0), # Remove batch dimension if only one sequence
        'attention_mask': encoding['attention_mask'].squeeze(0)
    }



def predict_with_sliding_window(
    long_sequence: str,
    model: AutoModelForSequenceClassification,
    tokenizer: AutoTokenizer,
    max_len: int = 400,
    stride: int = 1, # How many nucleotides to slide by each step
    dnabert_k: int = None,
    device: torch.device = torch.device("cpu")
) -> list:
    """
    Performs predictions on a long DNA sequence using a sliding window.

    Args:
        long_sequence (str): The full DNA sequence to scan.
        model (...): The loaded SpliceBERT model.
        tokenizer (...): The loaded tokenizer.
        max_len (int): The fixed input length for the model (e.g., 400).
        stride (int): How many nucleotides to slide the window by each step.
        dnabert_k (int, optional): K-mer size for DNA-BERT specific tokenization.
        device (torch.device): The device to run inference on.

    Returns:
        list: A list of tuples, where each tuple is (window_start_position, predicted_probability).
    """
    model.eval()
    predictions = []

    # Iterate through the sequence with the sliding window
    for start_pos in range(0, len(long_sequence) - max_len + 1, stride):
        end_pos = start_pos + max_len
        current_window_seq = long_sequence[start_pos:end_pos]

        # Prepare input for the current window
        inputs = prepare_single_sequence_input(
            current_window_seq,
            tokenizer,
            max_len=max_len,
            dnabert_k=dnabert_k
        )

        # Add batch dimension and move to device
        input_ids = inputs['input_ids'].unsqueeze(0).to(device)
        attention_mask = inputs['attention_mask'].unsqueeze(0).to(device)

        # Make prediction
        with torch.no_grad():
            with torch.amp.autocast(device.type): # Use the updated autocast syntax
                logits = model(input_ids=input_ids, attention_mask=attention_mask).logits.squeeze(1)
                probability = torch.sigmoid(logits).item()
        
        candidate_site_pos = start_pos + (max_len // 2)
        predictions.append((candidate_site_pos, logits.item(), probability))

    return predictions