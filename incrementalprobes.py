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
parser.add_argument("-b", "--batch", required=True, type=int, help="number of probes added per batch")
parser.add_argument("-m", "--max", required=True, type=int, help="maximum number of probes in panel")
parser.add_argument("-c", "--coverage", required=True, type=float, help="minimum desired coverage for 10th percentile of targets")
args=parser.parse_args()
#-------------------------------------------------------------------------------------

print("")

cmd='cp '+args.targets+' '+args.output+'_round_0_low_cov.fa'
subprocess.call(cmd,shell=True)

probes=0
cmd='touch '+args.output+'_final_probes.fa'
subprocess.call(cmd,shell=True)

r=1
loop=True

while loop==True:
    if r==1:
        cmd='cp '+args.targets+' '+args.output+'_round_0_low_cov.fa'
        subprocess.call(cmd,shell=True)
    else:
        cmd='python getlowcov.py -i '+args.output+'_round_'+str(r-1)+'_capture.pt -o '+args.output+'_round_'+str(r-1)
        subprocess.call(cmd,shell=True)
        
    cmd='python makeprobes.py -t '+args.output+'_round_'+str(r-1)+'_low_cov.fa -n '+str(args.batch)+' -o '+args.output+'_round_'+str(r)
    subprocess.call(cmd,shell=True)
    cmd='rm '+args.output+'_round_'+str(r-1)+'_low_cov.fa'
    subprocess.call(cmd,shell=True)
    
    with open(args.output+'_round_'+str(r)+'_probes.fa') as probesFile:
        probesFileLines=probesFile.readlines()
    newProbes=len([line for line in probesFileLines if line[0]=='>'])
    if newProbes>0:
        if r==1:
            cmd='python capture.py -t '+args.targets+' -p '+args.output+'_round_'+str(r)+'_probes.fa -o '+args.output+'_round_'+str(r)
            subprocess.call(cmd,shell=True)
        else:
            cmd='python capture.py -t '+args.targets+' -p '+args.output+'_round_'+str(r)+'_probes.fa -r '+args.output+'_round_'+str(r-1)+'_capture.pt -o '+args.output+'_round_'+str(r)
            subprocess.call(cmd,shell=True)
            cmd='rm '+args.output+'_round_'+str(r-1)+'_capture.pt'
            subprocess.call(cmd,shell=True)

        cmd='python stats.py -i '+args.output+'_round_'+str(r)+'_capture.pt -o '+args.output+'_round_'+str(r)
        subprocess.call(cmd,shell=True)
        cmd='rm '+args.output+'_round_'+str(r)+'_long.tsv'
        subprocess.call(cmd,shell=True)
        with open(args.output+'_round_'+str(r)+'_summary.tsv') as summaryFile:
            summaryFileLines=summaryFile.readlines()
        tenthPercentileCoverage=float(summaryFileLines[-1].split('\t')[4])

        with open(args.output+'_final_probes.fa') as probesFile:
            probesFileLines=probesFile.readlines()
        probes=len([line for line in probesFileLines if line[0]=='>'])
        if probes+newProbes>args.max:
            cmd='head -n '+str((args.max-probes)*2)+' '+args.output+'_round_'+str(r)+'_probes.fa >> '+args.output+'_final_probes.fa'
            subprocess.call(cmd,shell=True)
        else:
            cmd='cat '+args.output+'_round_'+str(r)+'_probes.fa >> '+args.output+'_final_probes.fa'
            subprocess.call(cmd,shell=True)
        cmd='rm '+args.output+'_round_'+str(r)+'_probes.fa'
        subprocess.call(cmd,shell=True)

        with open(args.output+'_final_probes.fa') as probesFile:
            probesFileLines=probesFile.readlines()
        probes=len([line for line in probesFileLines if line[0]=='>'])
        print('*'*50)
        print('END OF ROUND '+str(r))
        print('Probes added:',newProbes)
        print('Total probes:',probes)
        print('Tenth percentile of coverage:',tenthPercentileCoverage)
        print('*'*50)

        if tenthPercentileCoverage>=args.coverage:
            print('Desired coverage target achieved.')
            loop=False
        elif probes>=args.max:
            print('Maximum panel size reached.')
            loop=False
        else:
            r+=1

    else:
        print('No more probes can be designed.')
        loop=False

cmd='python capture.py -t '+args.targets+' -p '+args.output+'_final_probes.fa -o '+args.output+'_final'
subprocess.call(cmd,shell=True)
cmd='python stats.py -i '+args.output+'_final_capture.pt -o '+args.output+'_final'
subprocess.call(cmd,shell=True)
cmd='rm '+args.output+'_final_capture.pt'
subprocess.call(cmd,shell=True)

print()
print('Done.')
