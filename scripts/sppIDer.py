__author__ = 'Quinn'
__modified_by__ = 'Junhao Chen'

import argparse, multiprocessing, sys, re, subprocess, time
import os
from version import __version__

################################################################
# This script runs the full sppIDer pipeline.
# Before running this you must make your combination reference genome with the script combineRefGenomes.py
# This script will map short-read data to a combination reference genome and parse the outputs to create a summary of
# where and how well the reads map to each species in the combination reference genome.
#
# Program Requirements: bwa, samtools, bedtools, R, Rpackage ggplot2, Rpackage dplyr
# Input: Output name, Combination reference genome, fastq short-read files
#
################################################################

parser = argparse.ArgumentParser(description="Run full sppIDer")
parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--ref', help="Reference Genome, required", required=True)
parser.add_argument('--r1', help="Read1, required", required=True)
parser.add_argument('--r2', help="Read2, optional")
parser.add_argument('--byBP', help="Calculate coverage by basepair, optional, DEFAULT, can't be used with --byGroup", dest='bed', action='store_true')
parser.add_argument('--byGroup', help="Calculate coverage by chunks of same coverage, optional, can't be used with --byBP", dest='bed', action='store_false')
parser.add_argument('--seq-name', help="Output sequence name, optional")
parser.add_argument('--mq-threshold', type=int, help="Set MQ threshold, optional", default=3)
parser.add_argument('--cores', type=int, help="Set number of cores used in analysis, optional", default=(multiprocessing.cpu_count()//2))
parser.add_argument('--seq-type', choices=['PacBio', 'ONT'], help="Set sequence type (PacBio, ONT), optional")
parser.add_argument('--mapping-tool', choices=['bwa', 'minimap2'], help="Set mapping tool (default is bwa, optional is minimap2)", default='bwa')
parser.add_argument('--extract-species', help="Extract reads for specified species (comma-separated, e.g., Scer,Sbay)")
parser.add_argument('--extract-mq', type=int, default=30, help="Minimum MQ threshold for species extraction (default: 30)")
parser.add_argument('--extract-format', choices=['fastq', 'fastq.gz', 'list', 'both'], default='fastq.gz', help="Output format for extracted reads: fastq, fastq.gz (compressed), list, or both (default: fastq.gz)")
parser.add_argument('--skip-plot', action='store_true', help="Skip all plotting steps (MQ plot, depth plots) to save time")
parser.add_argument('--skip-depth', action='store_true', help="Skip depth calculation and plotting (bedtools and R depth analysis)")
parser.add_argument('--keep-sam', action='store_true', help="Keep SAM file after processing (default: delete to save space)")
parser.add_argument('--keep-intermediate', action='store_true', help="Keep intermediate BAM files (.view.bam)")
parser.set_defaults(bed=True)
args = parser.parse_args()

# docker vars
scriptDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '')
workingDir = os.path.join(os.getcwd(), '')

numCores = str(int(args.cores))
outputPrefix = args.out
refGen = args.ref
read1Name = args.r1
read2Name = args.r2 if args.r2 else None
start = time.time()

def calcElapsedTime(endTime):
    trackedTime = str()
    if 60 < endTime < 3600:
        min = int(endTime) // 60
        sec = int(endTime % 60)
        trackedTime = f"{min} mins {sec} secs"
    elif 3600 < endTime < 86400:
        hr = int(endTime) // 3600
        min = int((endTime % 3600) // 60)
        sec = int(endTime % 60)
        trackedTime = f"{hr} hrs {min} mins {sec} secs"
    elif 86400 < endTime < 604800:
        day = int(endTime) // 86400
        hr = int((endTime % 86400) // 3600)
        min = int((endTime % 3600) // 60)
        sec = int(endTime % 60)
        trackedTime = f"{day} days {hr} hrs {min} mins {sec} secs"
    elif 604800 < endTime:
        week = int(endTime) // 604800
        day = int((endTime % 604800) // 86400)
        hr = int((endTime % 86400) // 3600)
        min = int((endTime % 3600) // 60)
        sec = int(endTime % 60)
        trackedTime = f"{week} weeks {day} days {hr} hrs {min} mins {sec} secs"
    else:
        trackedTime = f"{int(endTime)} secs"
    return trackedTime

trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'w')
trackerOut.write(f"outputPrefix={args.out}\n")
trackerOut.write(f"ref={refGen}\n")
trackerOut.write(f"read1={read1Name}\n")
if read2Name: trackerOut.write(f"read2={read2Name}\n")
if not args.bed:
    trackerOut.write("coverage analysis option = by coverage groups, bedgraph format -bga\n")
else: trackerOut.write("coverage analysis option = by each base pair -d\n")
trackerOut.close()

########################## BWA ###########################
bwaOutName = outputPrefix + ".sam"
bwaOutFile = open(os.path.join(workingDir, bwaOutName), 'w')

match args.mapping_tool:
    case 'bwa':
        if args.seq_type and args.seq_type.lower() == 'pacbio':
            command = [args.mapping_tool, "mem", "-t", numCores, "-x", "pacbio", refGen, read1Name]
            print("Executing command:", " ".join(command))
            subprocess.call(command, stdout=bwaOutFile, cwd=workingDir)
        elif args.seq_type and args.seq_type.lower() == 'ont':
            command = [args.mapping_tool, "mem", "-t", numCores, "-x", "ont2d", refGen, read1Name]
            print("Executing command:", " ".join(command))
            subprocess.call(command, stdout=bwaOutFile, cwd=workingDir)
        else:
            command = [args.mapping_tool, "mem", "-t", numCores, refGen, read1Name]
            print("Executing command:", " ".join(command))
            subprocess.call(command, stdout=bwaOutFile, cwd=workingDir)
    case 'minimap2':
        if args.seq_type and args.seq_type.lower() == 'pacbio':
            command = ["minimap2", "-x", "map-pb", "-a", "-t", numCores, refGen, read1Name]
            print("Executing command:", " ".join(command))
            subprocess.call(command, stdout=bwaOutFile, cwd=workingDir)
        elif args.seq_type and args.seq_type.lower() == 'ont':
            command = ["minimap2", "-x", "map-ont", "-a", "-t", numCores, refGen, read1Name]
            print("Executing command:", " ".join(command))
            subprocess.call(command, stdout=bwaOutFile, cwd=workingDir)
        else:
            command = ["minimap2", "-a", "-t", numCores, refGen, read1Name]
            print("Executing command:", " ".join(command))
            subprocess.call(command, stdout=bwaOutFile, cwd=workingDir)
    case _:
        raise ValueError(f"Unsupported mapping tool: {args.mapping_tool}")

bwaOutFile.close()
print("BWA complete")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"BWA complete\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## samtools ###########################
bamSortOut = outputPrefix + ".sort.bam"

if not args.skip_depth:
    # Use pipe to avoid intermediate .view.bam file: sam -> view -> sort -> bam
    # samtools view filters by MQ, sort orders by position for bedtools
    print("Converting SAM to sorted BAM (using pipe)...")
    with open(os.path.join(workingDir, bamSortOut), 'w') as bamOut:
        view_proc = subprocess.Popen(
            ["samtools", "view", "-@", numCores, "-q", str(args.mq_threshold), "-bhSu", bwaOutName],
            stdout=subprocess.PIPE, cwd=workingDir
        )
        sort_proc = subprocess.Popen(
            ["samtools", "sort", "-@", numCores, "-o", bamSortOut, "-"],
            stdin=view_proc.stdout, cwd=workingDir
        )
        view_proc.stdout.close()
        sort_proc.communicate()
    print("SAMTOOLS complete (SAM -> sorted BAM via pipe)")
else:
    print("Skipped BAM conversion (--skip-depth)")

currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nSAMTOOLS complete\nElapsed time: {elapsedTime}")
trackerOut.close()

########################## parse SAM file ###########################
subprocess.call(["python3", os.path.join(scriptDir, "parseSamFile.py"), outputPrefix], cwd=workingDir)
print("Parsed SAM file")
currentTime = time.time() - start
elapsedTime = calcElapsedTime(currentTime)
print(f"Elapsed time: {elapsedTime}")
trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
trackerOut.write(f"\nParsed SAM\nElapsed time: {elapsedTime}")
trackerOut.close()

# Clean up SAM file to save disk space (default behavior)
samFilePath = os.path.join(workingDir, bwaOutName)
if not args.keep_sam and os.path.exists(samFilePath):
    os.remove(samFilePath)
    print(f"Removed {bwaOutName} to save disk space (use --keep-sam to retain)")

########################## plot MQ scores ###########################
if not args.skip_plot:
    subprocess.call(["Rscript", os.path.join(scriptDir, "MQscores_sumPlot.R"), outputPrefix], cwd=workingDir)
    print("Plotted MQ scores")
    currentTime = time.time() - start
    elapsedTime = calcElapsedTime(currentTime)
    print(f"Elapsed time: {elapsedTime}")
    trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
    trackerOut.write(f"\nMQ scores plotted\nElapsed time: {elapsedTime}")
    trackerOut.close()
else:
    print("Skipped MQ scores plotting (--skip-plot)")

########################## bedgraph Coverage ###########################
if not args.skip_depth:
    sortOut = bamSortOut
    if args.bed:
        bedOutD = outputPrefix + "-d.bedgraph"
        bedFileD = open(os.path.join(workingDir, bedOutD), 'w')
        subprocess.call(["genomeCoverageBed", "-d", "-ibam", sortOut], stdout=bedFileD, cwd=workingDir)
        bedFileD.close()
    else:
        bedOut = outputPrefix + ".bedgraph"
        bedFile = open(os.path.join(workingDir, bedOut), 'w')
        subprocess.call(["genomeCoverageBed", "-bga", "-ibam", sortOut], stdout=bedFile, cwd=workingDir)
        bedFile.close()
    print("bedgraph complete")
    currentTime = time.time() - start
    elapsedTime = calcElapsedTime(currentTime)
    print(f"Elapsed time: {elapsedTime}")
    trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
    trackerOut.write(f"\nbedgraph complete\nElapsed time: {elapsedTime}")
    trackerOut.close()

    ########################## average Bed ###########################
    if args.bed:
        subprocess.call(["Rscript", os.path.join(scriptDir, "meanDepth_sppIDer-d.R"), outputPrefix], cwd=workingDir)
    else:
        subprocess.call(["Rscript", os.path.join(scriptDir, "meanDepth_sppIDer-bga.R"), outputPrefix], cwd=workingDir)
    print("Found mean depth")
    currentTime = time.time() - start
    elapsedTime = calcElapsedTime(currentTime)
    print(f"Elapsed time: {elapsedTime}")
    trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
    trackerOut.write(f"\nFound mean depth\nElapsed time: {elapsedTime}")
    trackerOut.close()

    ########################## make plot ###########################
    if not args.skip_plot:
        subprocess.call(["Rscript", os.path.join(scriptDir, "sppIDer_depthPlot_forSpc.R"), outputPrefix], cwd=workingDir)
        subprocess.call(["Rscript", os.path.join(scriptDir, "sppIDer_depthPlot.R"), outputPrefix], cwd=workingDir)
        print("Plot complete")
        currentTime = time.time() - start
        elapsedTime = calcElapsedTime(currentTime)
        print(f"Elapsed time: {elapsedTime}")
        trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
        trackerOut.write(f"\nPlot complete\nElapsed time: {elapsedTime}\n")
        trackerOut.close()
    else:
        print("Skipped depth plotting (--skip-plot)")
else:
    print("Skipped depth calculation and plotting (--skip-depth)")

########################## extract reads by species (optional) ###########################
if args.extract_species:
    print(f"\nExtracting reads for species: {args.extract_species}")
    # Use absolute paths for fastq files
    fastq1_abs = os.path.abspath(read1Name)
    extractCmd = [
        "python3", os.path.join(scriptDir, "extractReadsBySpecies.py"),
        "--mq-file", os.path.join(workingDir, outputPrefix + "_MQ.txt"),
        "--fastq", fastq1_abs,
        "--species", args.extract_species,
        "--mq-min", str(args.extract_mq),
        "--output-format", args.extract_format
    ]
    if read2Name:
        fastq2_abs = os.path.abspath(read2Name)
        extractCmd.extend(["--fastq2", fastq2_abs])

    print("Executing command:", " ".join(extractCmd))
    result = subprocess.call(extractCmd, cwd=workingDir)

    if result == 0:
        print("Read extraction complete")
    else:
        print("Read extraction completed with errors")
    currentTime = time.time() - start
    elapsedTime = calcElapsedTime(currentTime)
    print(f"Elapsed time: {elapsedTime}")
    trackerOut = open(os.path.join(workingDir, outputPrefix + "_sppIDerRun.info"), 'a')
    trackerOut.write(f"\nRead extraction complete (species: {args.extract_species}, MQ>={args.extract_mq})\nElapsed time: {elapsedTime}\n")
    trackerOut.close()
