__author__ = 'Quinn'
__modified_by__ = 'Junhao Chen'

import sys, re, time
from pathlib import Path

################################################################
# This script will take the sam formatted output of bwa and parse for mapping quality and which genomes reads map to.
#
# Input: sam file of reads mapped to a combination reference genome
################################################################

# 获取脚本所在目录和当前工作目录
scriptDir = Path(__file__).resolve().parent
workingDir = Path.cwd()

inputName = sys.argv[1]
samName = inputName + ".sam"
outputName = inputName + "_MQ.txt"
outputLenName = inputName + "_chrLens.txt"
start = time.time()

speciesDict = {}
speciesDict["*"] = {}
speciesDict["*"][0] = {'count': 0, 'names': []}
speciesList = ['*']

# First pass: read header lines to get chromosome info, then process alignments
# Use line-by-line reading to avoid loading entire SAM file into memory
with open(workingDir / outputLenName, 'w', encoding='utf-8') as outputLen, \
     open(workingDir / samName, 'r', encoding='utf-8') as sam:

    for line in sam:
        if line.startswith('@SQ'):
            # Header line with sequence info
            headerInfo = line.strip().split('\t')
            chrInfo = headerInfo[1].split(":")[1]
            chrName = chrInfo.split("-")
            speciesName = chrName[0]
            chrNum = int(chrName[1])
            chrLen = headerInfo[2].split(":")[1]
            outputLen.write(chrInfo + "\t" + chrLen + "\n")
            if chrNum == 1:
                speciesList.append(speciesName)
                speciesDict[speciesName] = {}
                for i in range(0, 61):
                    speciesDict[speciesName][i] = {'count': 0, 'names': []}
        elif not line.startswith('@'):
            # Alignment line
            lineSplit = line.split('\t')
            if len(lineSplit) < 5:
                continue
            sequenceName = lineSplit[0]
            chrField = lineSplit[2]

            if chrField == '*':
                species = '*'
            else:
                species = chrField.split("-")[0]

            MQscore = int(lineSplit[4])

            if species not in speciesDict:
                speciesDict[species] = {}
            if MQscore not in speciesDict[species]:
                speciesDict[species][MQscore] = {'count': 0, 'names': []}

            speciesDict[species][MQscore]['count'] += 1
            speciesDict[species][MQscore]['names'].append(sequenceName)

with open(workingDir / outputName, 'w', encoding='utf-8') as output:
    output.write("Species\tMQscore\tcount\tSequenceNames\n")
    for species in speciesList:
        if species not in speciesDict:
            continue
        for score in sorted(speciesDict[species].keys()):
            count = speciesDict[species][score]['count']
            if count == 0:
                continue
            names = ",".join(speciesDict[species][score]['names'])
            output.write(f"{species}\t{score}\t{count}\t{names}\n")

currentTime = time.time() - start
print(f"{currentTime:.2f} secs")
