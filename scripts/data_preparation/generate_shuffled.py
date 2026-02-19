import random
from argparse import ArgumentParser
from typing import List

random.seed(42)

def append_record_to_container_for_fasta(sequence: str, identifier: str, container: List[str]):
    record = f">{identifier.rstrip()}\n{sequence.rstrip()}\n"
    container.append(record)

def shuffle_sequence(sequence: str) -> str:
    return "".join(random.sample(sequence, len(sequence)))

def main(input_path: str, output_path: str, num_shuffles: int):
    with open(input_path, "r") as f:
        fasta_lines = f.readlines()

    idx = 0
    new_records = []
    while idx < len(fasta_lines):
        header, sequence = fasta_lines[idx], fasta_lines[idx + 1]
        header = header.strip('>')
        # Append original sequence
        append_record_to_container_for_fasta(sequence, header, new_records)

        # Append shuffled sequences
        for i in range(num_shuffles):
            shuffled_seq = shuffle_sequence(sequence.strip())
            append_record_to_container_for_fasta(shuffled_seq, f"{header.rstrip()}_shuffled_{i+1}", new_records)

        idx += 2

    with open(output_path, "w") as fp:
        fp.writelines(new_records)

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Generate shuffled versions of sequences from a FASTA file."
    )
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input FASTA file containing sequences to be shuffled"
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to the output FASTA file where shuffled sequences will be written"
    )
    parser.add_argument(
        "-n", "--num-shuffles",
        type=int,
        default=10,
        help="Number of shuffled sequences to generate per original sequence (default: 10)"
    )
    args = parser.parse_args()

    main(args.input, args.output, args.num_shuffles)
    main(args.input, args.output, args.num_shuffles)
