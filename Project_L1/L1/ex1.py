'''
1. Make an app that is able to find the alphabet of a sequence of text. 
This sequence may be an RNA or an DNA or a protein sequence.
2. A DNA sequence is given: S= "ACGGGCATATGCGC." Make an app that is able
to show the percentage of the components from the alphabet of the sequence s. 
In other words, the input of the seq s and the output is the alphabet of the 
sequence and the percentage of each letter in the alphabet found in seq s.
3. Use the AI to adapt your current algorithm in order to make an app that
takes a FASTA file and read the sequence content from it and display the relative
percentages for the symbols present in the alphabet of sequences.
Note: FASTA represents a file format that contains DNA, RNA, or proteins sequence.
Thus, it contains the information for your input.
'''


def find_alphabet(sequence: str):
    return sorted(set(sequence))


def composition_percentage(sequence: str):
    length = len(sequence)
    percentages = {}
    for char in set(sequence):
        percentages[char] = (sequence.count(char) / length) * 100
    return percentages


def analyze_sequence(sequence: str):
    sequence = sequence.strip().upper()
    alphabet = find_alphabet(sequence)
    percentages = composition_percentage(sequence)

    print(f"\nSequence: {sequence}")
    print(f"Alphabet: {alphabet}")
    print("Composition (%):")
    for ch in alphabet:
        print(f"  {ch}: {percentages[ch]:.2f}%")


if __name__ == "__main__":
    S = "ACGGGCATATGCGC"
    print("Example with predefined DNA sequence:")
    analyze_sequence(S)