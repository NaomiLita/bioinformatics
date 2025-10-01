# def find_alphabet(sequence: str):
#     return sorted(set(sequence))
#
#
# def composition_percentage(sequence: str):
#     length = len(sequence)
#     percentages = {}
#     for char in set(sequence):
#         percentages[char] = (sequence.count(char) / length) * 100
#     return percentages
#
#
# def analyze_sequence(sequence: str):
#     sequence = sequence.strip().upper()
#     alphabet = find_alphabet(sequence)
#     percentages = composition_percentage(sequence)
#
#     print(f"\nSequence: {sequence}")
#     print(f"Alphabet: {alphabet}")
#     print("Composition (%):")
#     for ch in alphabet:
#         print(f"  {ch}: {percentages[ch]:.2f}%")
#
#
# if __name__ == "__main__":
#     S = "ACGGGCATATGCGC"
#     print("Example with predefined DNA sequence:")
#     analyze_sequence(S)

import os
from collections import Counter
import tkinter as tk
from tkinter import scrolledtext

fasta_file_path = "example.fasta"

def parse_fasta(path):
    sequences = []
    header = None
    seq_lines = []
    with open(path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    sequences.append((header, ''.join(seq_lines)))
                header = line[1:].strip() or "unknown"
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        sequences.append((header, ''.join(seq_lines)))
    return sequences


def analyze_sequence(sequence):
    seq = sequence.replace(" ", "").upper()
    length = len(seq)
    if length == 0:
        return None
    counts = Counter(seq)
    alphabet = sorted(counts.keys())
    percentages = {ch: (counts[ch] / length) * 100 for ch in alphabet}
    return length, alphabet, counts, percentages


def generate_output():
    output_text.delete('1.0', tk.END)

    if not os.path.exists(fasta_file_path):
        output_text.insert(tk.END, f"Error: file not found: {fasta_file_path}\n")
        return

    sequences = parse_fasta(fasta_file_path)
    for header, seq in sequences:
        res = analyze_sequence(seq)
        if res is None:
            output_text.insert(tk.END, f"> {header} -- (empty sequence)\n\n")
            continue
        length, alphabet, counts, percentages = res
        output_text.insert(tk.END, f"> {header}\n")
        output_text.insert(tk.END, f"Length: {length}\n")
        output_text.insert(tk.END, f"Alphabet: {alphabet}\n")
        output_text.insert(tk.END, "Composition:\n")
        for ch in alphabet:
            output_text.insert(tk.END, f"  {ch}: {counts[ch]} ({percentages[ch]:.2f}%)\n")
        output_text.insert(tk.END, "\n")

    if len(sequences) > 1:
        combined = "".join(seq for _, seq in sequences)
        res = analyze_sequence(combined)
        length, alphabet, counts, percentages = res
        output_text.insert(tk.END, "Combined analysis (all sequences):\n")
        output_text.insert(tk.END, f"Length: {length}\n")
        output_text.insert(tk.END, f"Alphabet: {alphabet}\n")
        output_text.insert(tk.END, "Composition:\n")
        for ch in alphabet:
            output_text.insert(tk.END, f"  {ch}: {counts[ch]} ({percentages[ch]:.2f}%)\n")
        output_text.insert(tk.END, "\n")

root = tk.Tk()
root.title("FASTA Sequence Analyzer")

generate_button = tk.Button(root, text="Generate Frequency Count", command=generate_output)
generate_button.pack(pady=10)
output_text = scrolledtext.ScrolledText(root, width=60, height=20)
output_text.pack(padx=10, pady=10)
root.mainloop()
