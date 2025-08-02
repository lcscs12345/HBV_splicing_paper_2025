import logging
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score
import numpy as np
import logging



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



def cal_metric_by_group(labels, preds, metric_fun, donors: list=None, acceptors: list=None, by_group: bool=True):
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
    
    
