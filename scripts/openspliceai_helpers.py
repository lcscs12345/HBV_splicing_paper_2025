import pandas as pd
import torch
from openspliceai.predict import predict


# description: Wrapper for OpenSpliceAI to generate and format splicing predictions from input sequences using PyTorch and pandas.


def openspliceai_predict(input_sequence, model, flanking_size):
    predictions = predict.predict(input_sequence, model, flanking_size)
    predictions_np = predictions.cpu().numpy() if isinstance(predictions, torch.Tensor) else predictions
    preds = pd.DataFrame(predictions_np)
    
    preds["position"] = preds.index -(flanking_size/2)
    preds["position"] = preds["position"].astype(int)
    length = len(input_sequence) - flanking_size
    preds = preds[(preds["position"]>=0) & (preds["position"]<=length)].copy()
    
    return preds