import torch
from transformers import AutoModel, AutoTokenizer

model_name = "facebook/esm2_t33_650M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name, dtype=torch.float16).eval()

sequence = "MKVLWAALLVTALAAGSLAEAAATA"  # Example AA sequence
inputs = tokenizer(sequence, return_tensors="pt", padding=True)

with torch.no_grad():
    outputs = model(**inputs)
    embeddings = outputs.last_hidden_state  # Shape: (1, seq_len+2, 1280) incl. BOS/EOS
