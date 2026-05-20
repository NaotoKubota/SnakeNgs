import os
import re

# Normalize the `samples` config entry.
#
# Accepts either:
#   - a list of sample names (legacy form), e.g. ["A", "B", "C"]
#   - a dict mapping a sample name to one or more run-level fastq prefixes
#     (multi-fastq form), e.g. {"sample1": ["XXX", "YYY"], "sample2": ["ZZZ"]}
#
# Returns a tuple (samples_list, samples_dict_or_None):
#   - samples_list: list of sample names used as the {sample} wildcard.
#   - samples_dict_or_None: the original dict when dict form is provided,
#     otherwise None (signals the legacy/no-merge path).
def normalize_samples(samples_cfg):
    if isinstance(samples_cfg, dict):
        samples_list = [str(x) for x in samples_cfg.keys()]
        samples_dict = {str(k): [str(v) for v in vs] for k, vs in samples_cfg.items()}
        return samples_list, samples_dict
    if isinstance(samples_cfg, (list, tuple)):
        return [str(x) for x in samples_cfg], None
    raise TypeError(
        f"config['samples'] must be a list or a dict, got {type(samples_cfg).__name__}"
    )

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
