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
parser.add_argument("-t", "--targets", type=str, required=True, help="path and name of the target input FASTA")
parser.add_argument("-o", "--output", type=str, required=True, help="job name to append to output files")
parser.add_argument("-k", "--size", nargs='?', default=120, type=int, help="length of Kmers (default=120)")
parser.add_argument("-s", "--step", nargs='?', default=1, type=int, help="Kmer step distance (default=1)")
parser.add_argument("-d", "--degen", nargs='?', default=0, type=int, help="Number of degenerate nucleotides permitted in a Kmer (default=0)")
parser.add_argument("-i", "--identity", nargs='?', default=0.90, type=float, help="Probe clustering identity (default=0.90 (90%))")
parser.add_argument("-n", "--number", nargs='?', default="max", type=str, help="Number of probe candidates returned (default=max)")
args=parser.parse_args()
#-------------------------------------------------------------------------------------

print("")

#---this section enumerates Kmers from the target input FASTA ------------------------
print("PROBETOOLS - MAKEPROBES [" + args.output + "]: Enumerating all " + '{:,}'.format(args.size) + "mers in " + args.targets + " ...")

targetSeqs={}
with open(args.targets) as targetFile:
    targetFileLines=targetFile.readlines()
for i in range(len(targetFileLines)):
    if targetFileLines[i][0]=='>':
        header=targetFileLines[i].lstrip('>').rstrip()
        if header in targetSeqs.keys():
            print("PROBETOOLS - MAKEPROBES [" + args.output + "]: WARNING --- Target " + header + " is duplicated in target dataset! Overwriting with latest sequence.")
        targetSeqs[header]=''
    else:
        targetSeqs[header]+=targetFileLines[i].rstrip()

outputFile=open(args.output + "_Kmers.fa", "w")
targets=0
kmers=0
for header,sequence in targetSeqs.items():
    if len(sequence)<args.size:
        print('PROBETOOLS - MAKEPROBES [' + args.output + ']: ERROR --- Target ' + header[:50] + ' is shorter than the desired probe length!')
    elif len(sequence)==0:
        print('PROBETOOLS - MAKEPROBES [' + args.output + ']: WARNING --- Target ' + header[:50] + ' has no sequence!')
    else:
        targets+=1
        for i in range(0, len(sequence)-args.size+1,args.step):
            k=sequence[i:i+args.size].upper()
            canonicalBases=sum([k.count(base) for base in list('ATGC')])
            #degenBases=sum([k.count(base) for base in  list('WSMKRYBVDHN')] 
            if canonicalBases>=args.size-args.degen:
                outputFile.write('>' + header[:50] + '_Kmer_' + str(i+1) + '\n')
                outputFile.write(k + '\n')
                kmers+=1            

print('PROBETOOLS - MAKEPROBES [' + args.output + ']: Enurmeration of '+'{:,}'.format(kmers)+' '+'{:,}'.format(args.size)+'-mers from '+'{:,}'.format(targets)+' target sequences in '+ args.targets +' complete.')
print()
outputFile.close()

#-------------------------------------------------------------------------------------

#---this sections call VSEARCH to clusters the Kmers ---------------------------------
cmd="vsearch --cluster_fast " + args.output + "_Kmers.fa --id " + str(args.identity) + " --centroids " + args.output + "_centroids.fa --fasta_width 0 --sizeout"
print("PROBETOOLS - MAKEPROBES [" + args.output + "]: Clustering " + args.output + " Kmers at " + str(args.identity*100) + "% identity...")
subprocess.call(cmd,shell=True)
cmd="rm " + args.output + "_Kmers.fa"
subprocess.call(cmd,shell=True)
print("PROBETOOLS - MAKEPROBES [" + args.output + "]: Kmer clustering finished.")
print()
#-------------------------------------------------------------------------------------

#---this section ranks the cluster centroids based on cluster size -------------------

clusterSizes={}
with open(args.output + "_centroids.fa", "r") as centroidsFile:
    centroidsFileLines=centroidsFile.readlines()
for i in range(len(centroidsFileLines)):
    if centroidsFileLines[i][0]=='>':
        clusterSize=int(centroidsFileLines[i].split('size=')[-1].rstrip())
    else:
        clusterSizes[centroidsFileLines[i].rstrip()]=clusterSize

print("PROBETOOLS - MAKEPROBES [" + args.output + "]: Ranking " + '{:,}'.format(len(clusterSizes)) + " probe candidates...")
cmd="rm " + args.output + "_centroids.fa"
subprocess.call(cmd,shell=True)

probeSeqs=[key for key,value in sorted(clusterSizes.items(), key=lambda item: item[1], reverse=True)]

#-------------------------------------------------------------------------------------

#---this section writes out the ranked probes to a FASTA file ------------------------
outputFile=open(args.output + '_probes.fa', "w")

if args.number=="max":
    maxProbes=len(probeSeqs)
else:
    maxProbes=int(args.number)

i=0
for probeSeq in probeSeqs[:maxProbes]:
    i+=1
    outputFile.write('>' + args.output + '_probe_' + str(i) + '\n')
    outputFile.write(probeSeq + '\n')

print("PROBETOOLS - MAKEPROBES [" + args.output + "]: " + '{:,}'.format(i) + " probes written to " + args.output + "_probes.fa")
print("PROBETOOLS - MAKEPROBES [" + args.output + "]: Finished.")
print("")

outputFile.close()
#-------------------------------------------------------------------------------------
