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
import numpy as np                                                                                                                                                                           
import pandas as pd

#---this section interprets the command line options when probetools.py is run----------                                                                                                                   
parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="path and name of the probetools capture results to be analyzed")
parser.add_argument("-o", "--output", type=str, required=True, help="job name to append to output files")
args=parser.parse_args()
#-------------------------------------------------------------------------------------                                                                                                                     

#--- this section opens the probetools results, calculates coverage and depth stats, and writes them to the output file--                                                                                  

print("PROBETOOLS - STATS [" + args.output + "]: Calculating statistics for file " + args.input + "...")

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

results=pd.DataFrame()
results['target']=list(targetSeqs.keys())
results['length']=results.apply(lambda row: len(targetSeqs[row.target]), axis=1)
results['Ns']=results.apply(lambda row: targetSeqs[row.target].count('N'), axis=1)
results['depth']=results.apply(lambda row: [depth for nucleotide,depth in zip(targetSeqs[row.target],targetDepths[row.target]) if nucleotide!='N'], axis=1)
results['0']=results.apply(lambda row: row['depth'].count(0), axis=1)
results['1']=results.apply(lambda row: row['depth'].count(1), axis=1)
results['2']=results.apply(lambda row: row['depth'].count(2), axis=1)
results['3']=results.apply(lambda row: row['depth'].count(3), axis=1)
results['4']=results.apply(lambda row: row['depth'].count(4), axis=1)
results['5+']=results.apply(lambda row: [depth>=5 for depth in row['depth']].count(True), axis=1)
cols=['1','2','3','4','5+']
results['covered']=results.apply(lambda row: row[cols].sum(), axis=1)
cols=['0','1','2','3','4','5+']
results['total']=results.apply(lambda row: row[cols].sum(), axis=1)
results['coverage']=results.apply(lambda row: round(100*row['covered']/row['total'],2), axis=1)
cols=['target','length','Ns','0','1','2','3','4','5+','covered','total','coverage']
results[cols].to_csv(args.output+'_long.tsv', sep='\t', index=False)

summaryFile=open(args.output + "_summary.tsv", "w")
summaryFile.write(args.output+'\n')
summaryFile.write("Coverage\t\t\t\t\t\t\t\tDepth\t\t\t\t\t\tTargets\t\n")
summaryFile.write("Mean\tStDev\tMin\t5%tile\t10%tile\tQ1\tMed\tQ3\tMax\t0\t1\t2\t3\t4\t5+\tFound\tTotal\t%\n")
summaryFile.write(str(round(results['coverage'].mean(),2))+'\t')
summaryFile.write(str(round(results['coverage'].std(),2))+'\t')
summaryFile.write(str(round(results['coverage'].min(),2))+'\t')
summaryFile.write(str(round(results['coverage'].quantile(0.05),2))+'\t')
summaryFile.write(str(round(results['coverage'].quantile(0.1),2))+'\t')
summaryFile.write(str(round(results['coverage'].quantile(0.25),2))+'\t')
summaryFile.write(str(round(results['coverage'].quantile(0.5),2))+'\t')
summaryFile.write(str(round(results['coverage'].quantile(0.75),2))+'\t')
summaryFile.write(str(round(results['coverage'].max(),2))+'\t')
summaryFile.write(str(round(100*results['0'].sum()/results['total'].sum(),2))+'\t')
summaryFile.write(str(round(100*results['1'].sum()/results['total'].sum(),2))+'\t')
summaryFile.write(str(round(100*results['2'].sum()/results['total'].sum(),2))+'\t')
summaryFile.write(str(round(100*results['3'].sum()/results['total'].sum(),2))+'\t')
summaryFile.write(str(round(100*results['4'].sum()/results['total'].sum(),2))+'\t')
summaryFile.write(str(round(100*results['5+'].sum()/results['total'].sum(),2))+'\t')
summaryFile.write(str(len(results[results['coverage']>0]))+'\t')
summaryFile.write(str(len(results))+'\t')
summaryFile.write(str(round(100*len(results[results['coverage']>0])/len(results),2))+'\n')
summaryFile.close()

print("PROBETOOLS - STATS [" + args.output + "]: Finished.")
print("")
