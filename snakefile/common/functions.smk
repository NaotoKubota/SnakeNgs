import os
import re

# Function to check if chrosome length is over 512 Mbp
def check_chromosome_length(star_index):
    chromosome_length_file = f"{star_index}/chrLength.txt" # This file should contain chromosome lengths per line in the format: length1
    if not os.path.exists(chromosome_length_file):
        raise FileNotFoundError(f"Chromosome length file not found: {chromosome_length_file}")
    with open(chromosome_length_file, 'r') as f:
        for line in f:
            length = int(line.strip())
            if length > 512 * 10**6:  # 512 Mbp
                return True
    return False
