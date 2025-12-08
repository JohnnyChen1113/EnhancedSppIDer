# EnhancedSppIDer

**Version: 0.2.1**

An enhanced version of [sppIDer](https://github.com/GLBRC/sppIDer) for species identification and read extraction from sequencing data.

## Features

- Support for Illumina, PacBio, and Oxford Nanopore sequencing data
- Flexible mapping tools: BWA or minimap2
- Species-specific read extraction with configurable quality thresholds
- Optimized pipeline with reduced memory usage and disk I/O
- Skip options for faster processing when plots are not needed

## Installation

### Dependencies

Install via conda (recommended):
```
conda install -c bioconda enhancedsppider
```

Or you can also install it step by step like:

```bash
conda create -n sppider python=3.10
conda activate sppider

# Core dependencies
conda install -c bioconda bwa samtools bedtools seqtk
conda install -c conda-forge biopython

# Optional: minimap2 for long-read mapping
conda install -c bioconda minimap2

# R dependencies for plotting
conda install -c conda-forge r-base r-ggplot2 r-dplyr

# Optional: kingfisher for downloading SRA data
conda install kingfisher
```

### Clone the repository

```bash
git clone https://github.com/JohnnyChen1113/EnhancedSppIDer.git
cd EnhancedSppIDer
```

## Quick Start

### Step 1: Prepare Reference Genomes

Create a key file (tab-separated) listing species names and their FASTA files. Example files (`Seub.fasta`, `Suva.fasta`, `SbayKey.txt`) are provided in the `examples/` directory:

```bash
# Example key file (You can find it in examples/SbayKey.txt)
Seub	Seub.fasta
Suva	Suva.fasta
```

Combine reference genomes:

```bash
cd examples
python scripts/combineRefGenomes.py --key SbayKey.txt --out Sbay
```

This generates:
- `Sbay` - Combined reference genome
- `Sbay.amb`, `.ann`, `.bwt`, `.pac`, `.sa` - BWA index files
- `Sbay.fai` - Samtools index
- `comboLength_Sbay.txt` - Contig lengths

### Step 2: Download Test Data (For Example)

You can download sequencing data using any tool of your choice. Here is an example using [kingfisher](https://github.com/wwood/kingfisher-download):

```bash
# Example: Download data from SRA using kingfisher
kingfisher get -r ERR1544719 -m aws-http -f fastq.gz --download-threads 8
```


### Step 3: Run sppIDer Pipeline

#### Basic usage, same as original sppIDer (Illumina paired-end):

```bash
python ../scripts/sppIDer.py \
    --out Sbay_out \
    --ref Sbay \
    --r1 ERR1544719_1.fastq.gz \
    --r2 ERR1544719_2.fastq.gz \
    --cores 8
```

#### With species extraction (in the case that you want get the sequance of each species beyond plots):

```bash
python ../scripts/sppIDer.py \
    --out Sbay_out \
    --ref Sbay \
    --r1 ERR1544719_1.fastq.gz \
    --r2 ERR1544719_2.fastq.gz \
    --cores 8 \
    --extract-species Seub,Suva \
    --extract-mq 30 \
    --extract-format fastq.gz
```

#### Fast mode (skip plotting):

```bash
python ../scripts/sppIDer.py \
    --out Sbay_out \
    --ref Sbay \
    --r1 ERR1544719_1.fastq.gz \
    --r2 ERR1544719_2.fastq.gz \
    --cores 8 \
    --skip-plot \
    --extract-species Seub,Suva
```

#### Long-read data (PacBio/ONT):

```bash
# PacBio with minimap2
python ../scripts/sppIDer.py \
    --out pacbio_run \
    --ref Sbay \
    --r1 pacbio_reads.fastq.gz \
    --seq-type PacBio \
    --mapping-tool minimap2 \
    --cores 12

# Oxford Nanopore
python ../scripts/sppIDer.py \
    --out ont_run \
    --ref Sbay \
    --r1 nanopore_reads.fastq.gz \
    --seq-type ONT \
    --mapping-tool minimap2 \
    --cores 12
```

### Step 4: Extract Reads by Species (Standalone)

**Note:** The original sppIDer does not support extracting read IDs or species-specific reads. This is a new feature in EnhancedSppIDer.

If you already have sppIDer output and want to extract reads separately:

```bash
# Extract reads for specific species from MQ file
python ../scripts/extractReadsBySpecies.py \
    --mq-file Sbay_out_MQ.txt \
    --fastq ERR1544719_1.fastq.gz \
    --fastq2 ERR1544719_2.fastq.gz \
    --species Seub,Suva \
    --mq-min 30 \
    --output-format fastq.gz

# Output files:
#   Seub_1.fastq.gz, Seub_2.fastq.gz
#   Suva_1.fastq.gz, Suva_2.fastq.gz
#   Seub_ids.txt, Suva_ids.txt (if --output-format list)
```

#### Extract from SAM file directly:

```bash
# Single-end reads
python ../scripts/extractReadsBySpecies.py \
    --sam-file Sbay_out.sam \
    --fastq ERR1544719_1.fastq.gz \
    --species Seub \
    --mq-min 30 \
    --output-format fastq.gz

# Paired-end reads with both FASTQ and ID list output
python ../scripts/extractReadsBySpecies.py \
    --sam-file Sbay_out.sam \
    --fastq ERR1544719_1.fastq.gz \
    --fastq2 ERR1544719_2.fastq.gz \
    --species Seub,Suva \
    --mq-min 30 \
    --output-format both

# Output files:
#   Seub_1.fastq.gz, Seub_2.fastq.gz, Seub_ids.txt
#   Suva_1.fastq.gz, Suva_2.fastq.gz, Suva_ids.txt
```

#### Output only ID lists (no FASTQ extraction):

```bash
python ../scripts/extractReadsBySpecies.py \
    --mq-file Sbay_out_MQ.txt \
    --species Seub,Suva \
    --mq-min 30 \
    --output-format list
```

## Command Reference

### sppIDer.py

| Option | Description | Default |
|--------|-------------|---------|
| `--out` | Output prefix (required) | - |
| `--ref` | Reference genome (required) | - |
| `--r1` | Read 1 FASTQ file (required) | - |
| `--r2` | Read 2 FASTQ file (paired-end) | - |
| `--cores` | Number of CPU cores | half of available |
| `--mq-threshold` | Mapping quality threshold | 3 |
| `--seq-type` | Sequence type: PacBio, ONT | - |
| `--mapping-tool` | Mapping tool: bwa, minimap2 | bwa |
| `--extract-species` | Species to extract (comma-separated) | - |
| `--extract-mq` | MQ threshold for extraction | 30 |
| `--extract-format` | Output format: fastq.gz, list, both | fastq.gz |
| `--skip-plot` | Skip all plotting steps | false |
| `--skip-depth` | Skip depth calculation and plotting | false |
| `--keep-sam` | Keep SAM file after processing | false |
| `--keep-intermediate` | Keep intermediate BAM files | false |
| `--byBP` | Calculate coverage by basepair (default) | true |
| `--byGroup` | Calculate coverage by groups | false |

### extractReadsBySpecies.py

| Option | Description | Default |
|--------|-------------|---------|
| `--mq-file` | Input _MQ.txt file from sppIDer | - |
| `--sam-file` | Input SAM file (alternative to --mq-file) | - |
| `--fastq` | Input FASTQ file | - |
| `--fastq2` | Input FASTQ file for read 2 | - |
| `--species` | Target species (comma-separated, required) | - |
| `--mq-min` | Minimum mapping quality | 0 |
| `--mq-max` | Maximum mapping quality | 60 |
| `--output-format` | Output: fastq.gz, list, both | fastq.gz |
| `--out` | Output directory | current directory |
| `--include-unmapped` | Include unmapped reads | false |
| `--no-seqtk` | Force Python extraction (slower) | false |

### combineRefGenomes.py

| Option | Description | Default |
|--------|-------------|---------|
| `--out` | Output filename (required) | - |
| `--key` | Key file with species and FASTA paths (required) | - |
| `--trim` | Minimum contig length to include | 0 |

## Output Files

### sppIDer.py outputs:

| File | Description |
|------|-------------|
| `{out}_MQ.txt` | Mapping quality scores by species |
| `{out}_chrLens.txt` | Chromosome lengths |
| `{out}.sort.bam` | Sorted BAM file |
| `{out}-d.bedgraph` | Coverage depth per base |
| `{out}_sppIDerRun.info` | Run log with timing |
| `{out}_MQsummary.pdf` | MQ score distribution plot |
| `{out}_depth.pdf` | Depth plot by species |
| `{species}_1.fastq.gz` | Extracted reads (if --extract-species) |

### extractReadsBySpecies.py outputs:

| File | Description |
|------|-------------|
| `{species}_1.fastq.gz` | Extracted R1 reads |
| `{species}_2.fastq.gz` | Extracted R2 reads (paired-end) |
| `{species}_ids.txt` | Read ID list |

## Performance Tips

1. **Use seqtk**: Install seqtk for 10x faster read extraction
2. **Skip plotting**: Use `--skip-plot` if you only need the data
3. **Skip depth**: Use `--skip-depth` for even faster runs (only MQ analysis)
4. **Adjust cores**: Use `--cores` to match your system

## Changelog

### v0.2.0
- Removed Docker dependency, works in conda environment
- Added species-specific read extraction (`--extract-species`)
- Added `--skip-plot` and `--skip-depth` options
- Optimized memory usage in SAM parsing (line-by-line reading)
- Added seqtk support for fast FASTQ extraction (with Python fallback)
- Auto-cleanup of intermediate files (SAM deleted by default)
- Support for compressed FASTQ output (.fastq.gz)
- Added `--version` flag to all main scripts
- Fixed CPU core type issue (float to int)
- Fixed file handle leaks in combineRefGenomes.py

### v0.1.0 (Original Enhanced Version)
- Added support for PacBio and ONT long-read sequencing
- Added minimap2 as alternative mapping tool
- Added `--cores` option for parallel processing

## License

MIT License

## Citation

If you use this tool, please cite the original sppIDer paper and this enhanced version.
