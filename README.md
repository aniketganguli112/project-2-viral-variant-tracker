# Viral Variant Tracker

A pipeline that takes SRA accessions of viral samples and tracks how
protein mutations change over time.

Designed for **SARS-CoV-2** (for example spike protein `S`). It reads
protein names directly from Nextclade output, so it works with any virus
that has a Nextclade dataset. **So far it has been tested on influenza A
(H1N1) HA samples**, chosen because they are small and quick to download.
SARS-CoV-2 runs are not yet validated (see Limitations).

## What it does

1. Downloads reads from SRA (`fasterq-dump`)
2. Quality control (`FastQC`) and trimming (`fastp`)
3. Alignment (`bwa mem`), sorting and indexing (`samtools`)
4. Coverage and low-coverage regions (`samtools depth`)
5. Variant calling and filtering (`bcftools`)
6. Consensus genome per sample (`bcftools consensus`)
7. Clade assignment and mutation calling (`Nextclade`)
8. Python analysis (`pandas`, `seaborn`):
   - QC filter (Nextclade status and coverage)
   - Mutation extraction for user-chosen proteins
   - Monthly mutation frequency table
   - Line plot and heatmap of frequency over time

## Repository structure

```
.
├── pipeline.sh            # reads -> consensus -> Nextclade
├── python_script/
│   └── script2.py         # frequency analysis and plots
├── example/
│   └── dates.csv          # example collection dates
├── refs/
│   └── h1n1pdm_ha_NC_026433.1.fasta   # test reference
├── environment.yml
├── .gitignore
└── README.md
```

## Requirements

Linux or WSL with conda.

```bash
conda env create -f environment.yml
conda activate project_2
```

Main tools: sra-tools, fastqc, fastp, bwa, samtools, bcftools, htslib,
nextclade, pandas, seaborn, matplotlib.

## Usage

### Step 1: run the pipeline

```bash
chmod +x pipeline.sh
./pipeline.sh
```

You will be asked for SRA accessions (separated by spaces), the path to
the reference FASTA, and the Nextclade dataset name.
The reference and the dataset **must match**.

| Virus | Reference | Nextclade dataset |
|---|---|---|
| SARS-CoV-2 | NC_045512.2 | `nextstrain/sars-cov-2/wuhan-hu-1/orfs` |
| Influenza A H1N1pdm HA (tested) | NC_026433.1 | `nextstrain/flu/h1n1pdm/ha` |

Outputs: `nc.tsv` (Nextclade results) and `all_consensus.fasta`.

### Step 2: run the analysis

Create `dates.csv` with columns `acc_no,dates` (`YYYY-MM-DD`).
See `example/dates.csv`.

```bash
python python_script/script2.py
```

You will be asked for:
- the root folder containing `nc.tsv` and `dates.csv`
- minimum coverage (0 to 1, e.g. `0.9`)
- protein names separated by space (`S` for SARS-CoV-2, `HA1 HA2` for influenza)

## Output

- Frequency table: fraction of samples carrying each mutation per month
- Line plot of mutation frequency over time
- Heatmap of mutations by month

## How frequency is calculated

```
frequency = samples with the mutation in a month / total samples in that month
```

The denominator counts all samples that passed QC, including samples
with no mutation in the chosen protein.

## Limitations

- Tested only on influenza A (H1N1) HA data. SARS-CoV-2 is untested.
- With few samples, frequencies can only be 0 or 1 and show no real trend.
  Use many samples across several months.
- Inputs are typed at prompts, so a sample sheet would be more robust.
- Nextclade column names differ between datasets, so some optional
  columns may be absent.

## Future work

- Sample-sheet input (accession, reference, dataset per sample)
- Automatic detection of rising mutations
- Full QC report
- Validation on SARS-CoV-2 data

## Author

Aniket Ganguli. Bioinformatics project.
