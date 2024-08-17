
'''
Snakefile for SUPPA2 diffSplice between Ref and Alt

Usage:
    snakemake -s SUPPA2_diffSplice.smk --configfile config.yaml --cores <int> --use-singularity
'''

import os
import pandas as pd
workdir: config['workdir']
absolute_path_to_workdir = os.path.abspath(config['workdir'])

'''experiment_table
sample	fastq	group
sample_01	/path/to/sample_01_1.fastq.gz,/path/to/sample_01_2.fastq.gz	Ref
sample_02	/path/to/sample_02_1.fastq.gz,/path/to/sample_02_2.fastq.gz	Ref
sample_03	/path/to/sample_03_1.fastq.gz,/path/to/sample_03_2.fastq.gz	Ref
sample_04	/path/to/sample_04_1.fastq.gz,/path/to/sample_04_2.fastq.gz	Alt
sample_05	/path/to/sample_05_1.fastq.gz,/path/to/sample_05_2.fastq.gz	Alt
sample_06	/path/to/sample_06_1.fastq.gz,/path/to/sample_06_2.fastq.gz	Alt
'''

# Read in the experiment table
def read_experiment_table(experiment_table):
    df = pd.read_csv(experiment_table, sep = "\t")
    samples = df['sample'].tolist()
    samples_Ref = df[df['group'] == 'Ref']['sample'].tolist()
    samples_Alt = df[df['group'] == 'Alt']['sample'].tolist()
    # Make dictionary with sample name as key and fastq file as value
    '''
    {
        sample1: {
            R1: /path/to/sample1_R1.fastq.gz,
            R2: /path/to/sample1_R2.fastq.gz
        },
        sample2: {
            R1: /path/to/sample2_R1.fastq.gz,
            R2: /path/to/sample2_R2.fastq.gz
        }
    }
    '''
    sample_fastq_dict = {sample: {'R1': os.path.join(absolute_path_to_workdir, fastq.split(',')[0]), 'R2': os.path.join(absolute_path_to_workdir, fastq.split(',')[1])} for sample, fastq in zip(df['sample'], df['fastq'])}
    return samples, samples_Ref, samples_Alt, sample_fastq_dict

samples, samples_Ref, samples_Alt, sample_fastq_dict = read_experiment_table(config['experiment_table'])
events = ["SE", "A3", "A5", "MX", "RI"]

rule all:
    input:
        quant_sf = expand("Salmon/{sample}/quant.sf", sample = samples),
        diff_splice_dpsi = expand("SUPPA2/results/diffSplice_{event}.dpsi", event = events),
        diff_splice_psivec = expand("SUPPA2/results/diffSplice_{event}.psivec", event = events)

rule salmon:
    wildcard_constraints:
        sample = "|".join(samples)
    container:
        "docker://combinelab/salmon:1.10.1"
    input:
        R1 = lambda wildcards: sample_fastq_dict[wildcards.sample]['R1'],
        R2 = lambda wildcards: sample_fastq_dict[wildcards.sample]['R2'],
        salmon_index = config['salmon_index']
    output:
        quant_sf = "Salmon/{sample}/quant.sf"
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/salmon_{sample}.txt"
    log:
        "log/salmon_{sample}.log"
    shell:
        """
        salmon quant \
        -i {input.salmon_index} \
        -l A \
        -1 {input.R1} \
        -2 {input.R2} \
        -p {threads} \
        -o Salmon/{wildcards.sample} \
        >& {log}
        """

rule merge_quant:
    input:
        quant_df_list = expand("Salmon/{sample}/quant.sf", sample=samples)
    output:
        merge_table = "SUPPA2/salmon_quant_merge.tsv"
    threads:
        1
    benchmark:
        "benchmark/merge_quant.txt"
    log:
        "log/merge_quant.log"
    run:
        import os
        import pandas as pd

        # Create an empty dataframe
        df = pd.DataFrame()

        # Read all the quant.sf files and concatenate them
        for f in input.quant_df_list:
            sample = os.path.basename(os.path.dirname(f))
            tmp = pd.read_csv(f, sep="\t", usecols=["Name", "TPM"])
            # Remove the version number from the gene ID
            tmp["Name"] = tmp["Name"].str.split(".").str[0]
            # Add the sample name to the dataframe
            tmp["sample"] = sample
            # Concatenate the dataframes
            df = pd.concat([df, tmp])

        # Pivot the dataframe
        df = df.pivot(index="Name", columns="sample", values="TPM")

        # Write a header with only the sample names
        header = "\t".join(df.columns)
        with open(output.merge_table, "w") as out:
            out.write(header + "\n")

        # Write the dataframe to the output file
        df.to_csv(output.merge_table, sep = "\t", index = True, header = False, mode="a")

rule generateEvents:
    container:
        "docker://naotokubota/suppa:2.3"
    output:
        suppa2_events_list = expand("SUPPA2/events/events_{event}_strict.ioe", event = events)
    threads:
        1
    benchmark:
        "benchmark/suppa2_generateEvents.txt"
    log:
        "log/suppa2_generateEvents.log"
    shell:
        """
        python /SUPPA-2.3/suppa.py generateEvents \
        -i {config[gtf]} \
        -o SUPPA2/events/events \
        -f ioe \
        -e SE SS MX RI
        """

rule psiPerEvent:
    wildcard_constraints:
        event = "|".join(events)
    container:
        "docker://naotokubota/suppa:2.3"
    input:
        merge_table = "SUPPA2/salmon_quant_merge.tsv",
        suppa2_events = "SUPPA2/events/events_{event}_strict.ioe"
    output:
        psi_table = "SUPPA2/results/results_{event}.psi"
    params:
        psi_output = "SUPPA2/results/results_{event}"
    threads:
        1
    benchmark:
        "benchmark/suppa2_psiPerEvent_{event}.txt"
    log:
        "log/suppa2_psiPerEvent_{event}.log"
    shell:
        """
        python /SUPPA-2.3/suppa.py psiPerEvent \
        -e {input.merge_table} \
        -i {input.suppa2_events} \
        -o {params.psi_output} \
        >& {log}
        """

rule separateGroupPsi:
    wildcard_constraints:
        event = "|".join(events)
    input:
        psi_table = "SUPPA2/results/results_{event}.psi"
    output:
        psi_Ref = "SUPPA2/results/psi_{event}_Ref.tsv",
        psi_Alt = "SUPPA2/results/psi_{event}_Alt.tsv"
    threads:
        1
    benchmark:
        "benchmark/suppa2_separateGroupPsi_{event}.txt"
    log:
        "log/suppa2_separateGroupPsi_{event}.log"
    run:
        import pandas as pd

        # Get header of the merge table
        header = ["id"] + list(pd.read_csv(input.psi_table, sep = "\t", nrows=1).columns)

        # Read the psi table
        psi_df = pd.read_csv(input.psi_table, sep = "\t", skiprows = 1, names = header)

        # Separate the psi table into two groups
        for group, samples in [("Ref", samples_Ref), ("Alt", samples_Alt)]:
            group_df = psi_df[["id"] + samples].fillna("nan")
            group_header = "\t".join(samples)
            with open(output[f"psi_{group}"], "w") as out:
                out.write(group_header + "\n")
            group_df.to_csv(output[f"psi_{group}"], sep = "\t", index = False, header = False, mode = "a")

rule separateGroupTpm:
    input:
        merge_table = "SUPPA2/salmon_quant_merge.tsv",
    output:
        tpm_Ref = "SUPPA2/results/tpm_Ref.tsv",
        tpm_Alt = "SUPPA2/results/tpm_Alt.tsv"
    threads:
        1
    benchmark:
        "benchmark/suppa2_separateGroupTpm.txt"
    log:
        "log/suppa2_separateGroupTpm.log"
    run:
        import pandas as pd

        # Get header of the merge table
        header = ["id"] + list(pd.read_csv(input.merge_table, sep = "\t", nrows=1).columns)

        # Read the merge table
        merge_df = pd.read_csv(input.merge_table, sep = "\t", skiprows = 1, names = header)

        # Separate the merge table into two groups
        for group, samples in [("Ref", samples_Ref), ("Alt", samples_Alt)]:
            group_df = merge_df[["id"] + samples].fillna("nan")
            group_header = "\t".join(samples)
            with open(output[f"tpm_{group}"], "w") as out:
                out.write(group_header + "\n")
            group_df.to_csv(output[f"tpm_{group}"], sep = "\t", index = False, header = False, mode = "a")

rule diffSplice:
    wildcard_constraints:
        event = "|".join(events)
    container:
        "docker://naotokubota/suppa:2.3"
    input:
        suppa2_events = "SUPPA2/events/events_{event}_strict.ioe",
        psi_Ref = "SUPPA2/results/psi_{event}_Ref.tsv",
        psi_Alt = "SUPPA2/results/psi_{event}_Alt.tsv",
        tpm_Ref = "SUPPA2/results/tpm_Ref.tsv",
        tpm_Alt = "SUPPA2/results/tpm_Alt.tsv"
    output:
        diff_splice_dpsi = "SUPPA2/results/diffSplice_{event}.dpsi",
        diff_splice_psivec = "SUPPA2/results/diffSplice_{event}.psivec"
    params:
        diff_splice_output = absolute_path_to_workdir + "/SUPPA2/results/diffSplice_{event}"
    threads:
        1
    benchmark:
        "benchmark/suppa2_diffSplice_{event}.txt"
    log:
        "log/suppa2_diffSplice_{event}.log"
    shell:
        """
        python /SUPPA-2.3/suppa.py diffSplice \
		--method empirical \
		--input {input.suppa2_events} \
		--psi \
		{input.psi_Ref} \
		{input.psi_Alt} \
		--tpm \
		{input.tpm_Ref} \
		{input.tpm_Alt} \
		-gc \
		-o {params.diff_splice_output} \
        --mode DEBUG \
        >& {log}
        """
