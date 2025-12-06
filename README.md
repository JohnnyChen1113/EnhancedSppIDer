# EnhancedSppIDer

**Version: 0.2.0**

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
pip install kingfisher
```

### Clone the repository

```bash
git clone https://github.com/your-repo/EnhancedSppIDer.git
cd EnhancedSppIDer
```

## Quick Start

### Step 1: Prepare Reference Genomes

Create a key file (tab-separated) listing species names and their FASTA files:

```bash
# Create key file
cat > SbayKey.txt << 'EOF'
Scer	S288c.fasta
Sbay	GCA_016858285.1.fasta
EOF
```

Combine reference genomes:

```bash
python scripts/combineRefGenomes.py --key SbayKey.txt --out Sbay_combo.fasta
```

This generates:
- `Sbay_combo.fasta` - Combined reference genome
- `Sbay_combo.fasta.amb`, `.ann`, `.bwt`, `.pac`, `.sa` - BWA index files
- `Sbay_combo.fasta.fai` - Samtools index
- `comboLength_Sbay_combo.fasta.txt` - Contig lengths

### Step 2: Download Test Data (Optional)

```bash
# Download example data from SRA
kingfisher get -r ERR1544719 -m aws-http -f fastq.gz --download-threads 8
```

### Step 3: Run sppIDer Pipeline

#### Basic usage (Illumina paired-end):

```bash
python scripts/sppIDer.py \
    --out test_run \
    --ref Sbay_combo.fasta \
    --r1 ERR1544719_1.fastq.gz \
    --r2 ERR1544719_2.fastq.gz \
    --cores 8
```

#### With species extraction:

```bash
python scripts/sppIDer.py \
    --out test_run \
    --ref Sbay_combo.fasta \
    --r1 ERR1544719_1.fastq.gz \
    --r2 ERR1544719_2.fastq.gz \
    --cores 8 \
    --extract-species Scer,Sbay \
    --extract-mq 30 \
    --extract-format fastq.gz
```

#### Fast mode (skip plotting):

```bash
python scripts/sppIDer.py \
    --out test_run \
    --ref Sbay_combo.fasta \
    --r1 ERR1544719_1.fastq.gz \
    --r2 ERR1544719_2.fastq.gz \
    --cores 8 \
    --skip-plot \
    --extract-species Scer,Sbay
```

#### Long-read data (PacBio/ONT):

```bash
# PacBio with minimap2
python scripts/sppIDer.py \
    --out pacbio_run \
    --ref Sbay_combo.fasta \
    --r1 pacbio_reads.fastq.gz \
    --seq-type PacBio \
    --mapping-tool minimap2 \
    --cores 12

# Oxford Nanopore
python scripts/sppIDer.py \
    --out ont_run \
    --ref Sbay_combo.fasta \
    --r1 nanopore_reads.fastq.gz \
    --seq-type ONT \
    --mapping-tool minimap2 \
    --cores 12
```

### Step 4: Extract Reads by Species (Standalone)

If you already have sppIDer output and want to extract reads separately:

```bash
# Extract reads for specific species from MQ file
python scripts/extractReadsBySpecies.py \
    --mq-file test_run_MQ.txt \
    --fastq ERR1544719_1.fastq.gz \
    --fastq2 ERR1544719_2.fastq.gz \
    --species Scer,Sbay \
    --mq-min 30 \
    --output-format fastq.gz

# Output files:
#   Scer_1.fastq.gz, Scer_2.fastq.gz
#   Sbay_1.fastq.gz, Sbay_2.fastq.gz
#   Scer_ids.txt, Sbay_ids.txt (if --output-format both)
```

#### Extract from SAM file directly:

```bash
python scripts/extractReadsBySpecies.py \
    --sam-file test_run.sam \
    --fastq reads.fastq.gz \
    --species Scer \
    --mq-min 30 \
    --output-format fastq.gz
```

#### Output only ID lists (no FASTQ extraction):

```bash
python scripts/extractReadsBySpecies.py \
    --mq-file test_run_MQ.txt \
    --species Scer,Sbay \
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
| `--extract-format` | Output format: fastq, fastq.gz, list, both | fastq.gz |
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
| `--output-format` | Output: fastq, fastq.gz, list, both | both |
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
