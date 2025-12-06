#!/usr/bin/env python
import sys, subprocess, argparse, os
from Bio import SeqIO
from version import __version__

################################################################
# This script will create the combination reference genome and all additional files needed to run the full sppIDer script.
#
# Program Requirements: bwa, samtools
# Input: output name prefix, text file key of genome names, and optional desired minimum length of contigs included (default is to keep all contigs)
# The text file key must be a tab seperated list of unique name for each genome and the actual fasta name for that genome, e.g.
# Saccharomyces_cerevisaie   S288c.fasta
# Saccharomyces_paradoxus    GCA_002079055.1.fasta
#
################################################################

# Use current working directory instead of docker path
workingDir = os.path.join(os.getcwd(), '')

parser = argparse.ArgumentParser(description="Combine desired reference genomes")
parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--key', help="Key to reference genomes, required", required=True)
parser.add_argument('--trim', type=int, help="Trim any contigs smaller than this value (for genomes in many small contigs)", default=0)
args = parser.parse_args()

comboGenomeName = args.out
keyFileName = args.key
trimLength = args.trim

# Check input files exist
keyFilePath = os.path.join(workingDir, keyFileName)
if not os.path.exists(keyFilePath):
    print(f"Error: Key file not found: {keyFilePath}", file=sys.stderr)
    sys.exit(1)

comboTotalLen = 0

with open(os.path.join(workingDir, "comboLength_" + comboGenomeName + ".txt"), 'w') as lengthFile, \
     open(os.path.join(workingDir, comboGenomeName), 'w') as outGenome, \
     open(keyFilePath, 'r') as keyFile:

    lengthFile.write(comboGenomeName + "\tTrimmed contigs >" + str(trimLength) + "\n")

    for line in keyFile:
        line = line.strip()
        if not line:
            continue

        parts = line.split('\t')
        if len(parts) < 2:
            print(f"Warning: Skipping invalid line: {line}", file=sys.stderr)
            continue

        uniID = parts[0]
        genomeName = parts[1]
        genomePath = os.path.join(workingDir, genomeName)

        if not os.path.exists(genomePath):
            print(f"Error: Genome file not found: {genomePath}", file=sys.stderr)
            sys.exit(1)

        sumGenomeLen = 0
        counter = 0

        with open(genomePath, 'r') as fasta:
            for seq_record in SeqIO.parse(fasta, "fasta"):
                if len(seq_record.seq) >= trimLength:
                    outGenome.write(">" + uniID + "-" + str(counter + 1) + "\n")
                    outGenome.write(str(seq_record.seq) + "\n")
                    lengthFile.write(uniID + "-" + str(counter + 1) + "\t" + str(len(seq_record.seq)) + "\n")
                    sumGenomeLen += len(seq_record.seq)
                    counter += 1

        comboTotalLen += sumGenomeLen

        # Format genome length for human readability
        if sumGenomeLen > 1000000000:
            strGenomeLen = f"{sumGenomeLen / 1e9:.2f} Gb"
        elif sumGenomeLen > 1000000:
            strGenomeLen = f"{sumGenomeLen / 1e6:.2f} Mb"
        elif sumGenomeLen > 1000:
            strGenomeLen = f"{sumGenomeLen / 1e3:.2f} Kb"
        else:
            strGenomeLen = f"{sumGenomeLen} bp"

        lengthFile.write(uniID + "-totalGenome\t" + strGenomeLen + "\n")
        print(f"  {uniID}: {counter} contigs, {strGenomeLen}")

    # Format total combo length
    if comboTotalLen > 1000000000:
        strComboLen = f"{comboTotalLen / 1e9:.2f} Gb"
    elif comboTotalLen > 1000000:
        strComboLen = f"{comboTotalLen / 1e6:.2f} Mb"
    elif comboTotalLen > 1000:
        strComboLen = f"{comboTotalLen / 1e3:.2f} Kb"
    else:
        strComboLen = f"{comboTotalLen} bp"

    lengthFile.write("Combo-totalGenome\t" + strComboLen + "\n")

print(f"Combined genome: {strComboLen}")
print("Building BWA index...")
result = subprocess.call(["bwa", "index", comboGenomeName], cwd=workingDir)
if result != 0:
    print("Error: BWA indexing failed", file=sys.stderr)
    sys.exit(1)

print("Building samtools index...")
result = subprocess.call(["samtools", "faidx", comboGenomeName], cwd=workingDir)
if result != 0:
    print("Error: samtools indexing failed", file=sys.stderr)
    sys.exit(1)

print("Done!")
