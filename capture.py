# ProbeTools                                                                                                                                                                                               
#                                                                                                                                                                                                          
# Kevin Kuchinski                                                                                                                                                                                          
# University of British Columbia, Department of Pathology and Laboratory Medicine                                                                                                                          
# British Columbia Centre for Disease Control, Public Health Laboratory, Prystajecky Lab                                                                                                                   
# kevin.kuchinski@bccdc.ca                                                                                                                                                                                 
#                                                                                                                                                                                                          
# Natalie Prystajecky, PhD SCCM (Env)                                                                                                                                                                      
# Program Head for Environmental Microbiology, BCCDC Public Health Laboratory                                                                                                                              
# Co-program Head for Molecular Microbiology and Genomics, BCCDC Public Health Laboratory                                                                                                                  
# Clinical Assistant Professor, Pathology & Laboratory Medicine, UBC                                                                                                                                      
# 655 West 12th Avenue                                                                                                                                                                                    
# Vancouver, BC, V5Z 4R4                                                                                                                                                                                  
# Canada                                                                                                                                                                                                   
# natalie.prystajecky@bccdc.ca 

import argparse
import subprocess

#---this section interprets the command line options when probetools.py is run -------
parser=argparse.ArgumentParser()
parser.add_argument("-p", "--probes", type=str, required=True, help="path and name of the probe FASTA")
parser.add_argument("-t", "--targets", type=str, required=True, help="path and name of the target FASTA")
parser.add_argument("-o", "--output", type=str, required=True, help="job name appended to output files")
parser.add_argument("-r", "--previous", nargs='?', default="", type=str, help="path and name of previous capture results for appending new results")
parser.add_argument("-i", "--identity", nargs='?', default=90, type=float, help="BLAST identity threshold (default=90 (90%))")
parser.add_argument("-l", "--length", nargs='?', default=60, type=int, help="BLAST alignment length threshold (default=60 (60 bp))")
args=parser.parse_args()

#---this section BLASTs the probes against the targets ------------------------------
print("")
print("PROBETOOLS - CAPTURE [" + args.output + "]: Aligning " + args.probes + " to " + args.targets + "...")
cmd="makeblastdb -in " + args.targets + " -dbtype nucl"
subprocess.call(cmd,shell=True)
cmd="blastn -db " + args.targets + " -query " + args.probes + " -num_threads 8 -max_target_seqs 100000 -outfmt 6 | awk '$3>=" + str(args.identity) + " {print}' | awk '$4>=" + str(args.length)
cmd+=" {print}' > " + args.output + "_BLAST.tsv"
subprocess.call(cmd,shell=True)
for suffix in ['.nhr','.nin','.nsq']:
    cmd="rm "+args.targets+suffix
    subprocess.call(cmd,shell=True)

print("PROBETOOLS - CAPTURE [" + args.output + "]: Finished aligning probes to targets.")

#---this section extracts the header names and sequences from the targets file and sets default probe depths to zero -------
print("PROBETOOLS - CAPTURE [" + args.output + "]: Extracting target names and sequences from " + args.targets + "...")
targetSeqs={}
with open(args.targets) as targetFile:
    targetFileLines=targetFile.readlines()
for i in range(len(targetFileLines)):
    if targetFileLines[i][0]=='>':
        header=targetFileLines[i].lstrip('>').rstrip()
        if header in targetSeqs.keys():
            print("PROBETOOLS - CAPTURE [" + args.output + "]: WARNING --- Target " + header + " is duplicated in target dataset! Overwriting with latest sequence.")
        targetSeqs[header]=''
    else:
        targetSeqs[header]+=targetFileLines[i].rstrip()

targetDepths={}
for header in targetSeqs.keys():
    targetDepths[header]=[0]*len(targetSeqs[header])

print("PROBETOOLS - CAPTURE [" + args.output + "]: Total targets extracted from " + args.targets + ": " + '{:,}'.format(len(targetSeqs)))

#---this section loads existing probe depth data -------
if args.previous!='':
    print("PROBETOOLS - CAPTURE [" + args.output + "]: Loading previous capture results from " + args.previous + "...")
    with open(args.previous) as previousFile:
        previousFileLines=previousFile.readlines()
    for i in range(len(previousFileLines)):
        if previousFileLines[i][0]=='>':
            header=previousFileLines[i].lstrip('>').rstrip()
            ignore=False
        elif previousFileLines[i][0]=='$':
            if header in targetSeqs.keys() and previousFileLines[i].lstrip('$').rstrip()!=targetSeqs[header]:
                print("PROBETOOLS - CAPTURE [" + job + "]: WARNING --- Sequence for target " + header + " in previous results does not match sequence in " + target + "! Ignoring previous results.")
                ignore=True
            elif header not in targetSeqs.keys():
                targetSeqs[header]=previousFileLines[i].lstrip('$').rstrip()
        elif previousFileLines[i][0]=='#':
            if ignore==False:
                targetDepths[header]=[int(depth) for depth in previousFileLines[i].lstrip('#').rstrip().split(',')]

#---this section parses the BLAST results and asigns probe depth to targets -------
print("PROBETOOLS - CAPTURE [" + args.output + "]: Extracting probe depth from alignment of probes against targets...")
BLASTresults=open(args.output + "_BLAST.tsv", "r")
currentLine=BLASTresults.readline()
while currentLine!="":
    currentTarget=currentLine.split("\t")[1].rstrip()
    alnStart=int(currentLine.split("\t")[8].rstrip())-1
    alnEnd=int(currentLine.split("\t")[9].rstrip())-1
    alnStartEnd=[alnStart,alnEnd]
    for i in range(min(alnStartEnd),max(alnStartEnd)+1):
        targetDepths[currentTarget][i]+=1
    currentLine=BLASTresults.readline()
BLASTresults.close()
cmd="rm " + args.output + "_BLAST.tsv"
subprocess.call(cmd,shell=True)

#---this section writes the target headers, sequences, and probe depths to the output file  -------
print("PROBETOOLS - CAPTURE [" + args.output + "]: Writing target headers, sequences, and probe depths to " + args.output+'_capture.pt' + " ...")
outputFile=open(args.output + '_capture.pt', "w")
for header in targetSeqs.keys():
    outputFile.write(">" + header + "\n")
    outputFile.write("$" + targetSeqs[header].upper() + "\n")
    outputFile.write("#" + ','.join(str(i) for i in targetDepths[header]) + "\n")

print("PROBETOOLS - CAPTURE [" + args.output + "]: Capture finished.")
print("")
