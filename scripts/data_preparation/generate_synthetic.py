import random
from collections import Counter
import argparse

random.seed(42)


def read_fasta(file_path):
    """Reads a FASTA file and returns a list of sequences."""
    sequences = []
    with open(file_path, 'r') as file:
        seq = ""
        for line in file:
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences


def get_aa_frequencies(sequences):
    """Computes amino acid frequencies across all sequences."""
    aa_counts = Counter()
    total_aa = 0
    for seq in sequences:
        aa_counts.update(seq)
        total_aa += len(seq)

    return {aa: count / total_aa for aa, count in aa_counts.items()}


def generate_synthetic_sequences(n, length_distribution, aa_frequencies, mode):
    """Generates N synthetic sequences based on the given length distribution and amino acid frequencies or uniform sampling."""
    synthetic_sequences = []
    amino_acids = list(aa_frequencies.keys())
    if mode == "random":
        probabilities = [1 / len(amino_acids)] * len(amino_acids)
    else:
        probabilities = list(aa_frequencies.values())

    for _ in range(n):
        seq_length = random.choice(length_distribution)
        synthetic_seq = ''.join(random.choices(amino_acids, probabilities, k=seq_length))
        synthetic_sequences.append(synthetic_seq)

    return synthetic_sequences


def write_fasta(sequences, output_file, mode):
    """Writes the generated sequences to a FASTA file."""
    id = 'synthetic_random' if mode == "random" else 'synthetic_realistic'
    with open(output_file, 'w') as file:
        for i, seq in enumerate(sequences):
            file.write(f">{id}_seq_{i + 1}\n")
            file.write(f"{seq}\n")


def main(input_fasta, output_fasta, n, mode):
    sequences = read_fasta(input_fasta)
    length_distribution = [len(seq) for seq in sequences]
    aa_frequencies = get_aa_frequencies(sequences)
    synthetic_sequences = generate_synthetic_sequences(n, length_distribution, aa_frequencies, mode)
    write_fasta(synthetic_sequences, output_fasta, mode)
    print(f"Generated {n} synthetic sequences and saved to {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic amino acid sequences.")
    parser.add_argument("input_fasta", type=str, help="Path to input FASTA file.")
    parser.add_argument("output_fasta", type=str, help="Path to output FASTA file.")
    parser.add_argument("-n", "--num-seqs",  type=int, help="Number of synthetic sequences to generate.")
    parser.add_argument("--mode", type=str, choices=["random", "realistic"], default="realistic",
                        help="Sampling mode: 'random' for equal amino acid probability, 'realistic' for frequency-based.")
    args = parser.parse_args()

    main(args.input_fasta, args.output_fasta, args.num_seqs, args.mode)
