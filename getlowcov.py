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

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="path and name of the probetools capture results to be analyzed")
parser.add_argument("-o", "--output", type=str, required=True, help="job name to append to output files")
parser.add_argument("-d", "--depth", nargs='?', default=0, type=int, help="lowcov will only return areas not exceeding this probe depth (default=0)")
parser.add_argument("-l", "--length", nargs='?', default=40, type=int, help="lowcov will only return areas of at least this length (default=40)")
parser.add_argument("-w", "--window", nargs='?', default=120, type=int, help="if the low coverage area is shorter than this value, lowcov will take a window of this size surrounding the lowcov area (default=120)")

args=parser.parse_args()

print("PROBETOOLS - GETLOWCOV [" + args.output + "]: Extracting sequences and probe coverage from file " + args.input + "...")

targetSeqs={}
targetDepths={}
with open(args.input, "r") as resultsFile:
    resultsLines=resultsFile.readlines()
for i in range(len(resultsLines)):
    if resultsLines[i][0]=='>':
        header=resultsLines[i].lstrip('>').rstrip()
    elif resultsLines[i][0]=='$':
        targetSeqs[header]=resultsLines[i].lstrip('$').rstrip()
    elif resultsLines[i][0]=='#':
        targetDepths[header]=[int(depth) for depth in resultsLines[i].lstrip('#').rstrip().split(',')]

outputFile=open(args.output+'_low_cov.fa', 'w')

for header in targetSeqs.keys():
    if len(targetSeqs[header])<args.window:
        print("PROBETOOLS - GETLOWCOV [" + args.output + "]: WARNING --- Window size exceeds length of target sequence " + header + " ! Ignoring target sequence.")
    else:
        lowCovSeqs=''.join([base if depth<=args.depth else 'X' for base,depth in zip(targetSeqs[header],targetDepths[header])]).split('X')
        lowCovSeqs=[seq for seq in lowCovSeqs if seq!='']
        lowCovSeqs=[seq for seq in lowCovSeqs if len(seq)>=args.length]
        shortLowCovSeqs=[seq for seq in lowCovSeqs if len(seq)<args.window]
        lowCovSeqs=[seq for seq in lowCovSeqs if seq not in shortLowCovSeqs]
        for seq in shortLowCovSeqs:
            seqStart=targetSeqs[header].index(seq)
            seqEnd=seqStart+len(seq)
            leftExtend=int((args.window-len(seq))/2)
            rightExtend=args.window-len(seq)-leftExtend
            extendedStart=seqStart-leftExtend
            extendedEnd=seqEnd+rightExtend
            if extendedStart<0:
                extendedEnd=extendedEnd-extendedStart
                extendedStart=extendedStart-extendedStart
            if extendedEnd>len(targetSeqs[header]):
                extendedStart=extendedStart-(extendedEnd-len(targetSeqs[header]))
                extendedEnd=extendedEnd-(extendedEnd-len(targetSeqs[header]))
            lowCovSeqs.append(targetSeqs[header][extendedStart:extendedEnd])
        i=0
        for seq in lowCovSeqs:
            outputFile.write('>'+header+'-'+str(i)+'\n')
            outputFile.write(seq+'\n')
            i+=1

outputFile.close()
print("PROBETOOLS - GETLOWCOV [" + args.output + "]: Finished.")
print("")
