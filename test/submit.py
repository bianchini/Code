#!/usr/bin/env python
import commands
import time
import re
import os
import string
from os import listdir
from os.path import isfile, join

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

path = '/pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/'
sample = 'TTJetsMG'
csvFile = 'TTbar_excl'
rnd_list = [0]
opt_list = [2]

################################################################################

files_per_job = 4

splits = []
dirnames = {
    'TTbar' : [path+'TT_TuneCUETP8M1_13TeV-powheg-pythia8/VHBB_HEPPY_V14_TT_TuneCUETP8M1_13TeV-powheg-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_212301/0000', 0],
    'TTJetsMG' : [path+'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_212453/0000/', 1],
    'WJets' : [path+'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_220559/0000/', 2],
    'ttH' : [path+'ttHTobb_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V14_ttHTobb_M125_13TeV_powheg_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_084144/0000/', 3],
}

dirname = dirnames[sample][0]

onlyfiles = [f for f in listdir(dirname) if isfile(join(dirname, f))]

files = []
for fl in onlyfiles:
    if 'tree' not in fl:
        continue
    fl1 = fl.lstrip('tree_') 
    fl2 = fl1.rstrip('.root') 
    files.append(int(fl2))
files_s = sorted(files)
fill = []
for k in range(len(files_s)):                                                      
    k = k+1
    if k % files_per_job!=0 and k in files_s and k+1 in files_s:
        fill.append(k)
    else:
        if k % files_per_job==0 and k in files_s:
            fill.append(k)
            splits.append(fill)   
            fill = []
                             
print splits
print "Total number of jobs: ", len(splits)





count = 0
for part in splits:
    count += 1
    for rnd in rnd_list:
        for option in opt_list:

            print "Job: ", count, ", rnd: ", rnd, ", option: ", option, ", files: ", part

            scriptName = 'job_'+str(rnd)+'_'+str(option)+'_p'+str(count)+'.sh'
            jobName    = 'job_'+str(rnd)+'_'+str(option)+'_p'+str(count)
            f = open(scriptName,'w')
            f.write('#!/bin/bash\n\n')
            f.write('cd /shome/bianchi/VHbb_3/CMSSW_7_4_15/src/TTH/MEIntegratorStandalone/test\n')
            f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
            f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
            f.write('eval `scramv1 runtime -sh`\n')
            f.write('\n')    
            f.write('../bin/btag_analyzer '+sample+' '+str(rnd)+' '+csvFile+' '+str(part[0])+' '+str(part[len(part)-1])+' '+str( dirnames[sample][1] )+' '+str(option)+'\n')
            f.close()
            os.system('chmod +x '+scriptName)
            
            submitToQueue = 'qsub -V -cwd -q all.q -N '+jobName+' '+scriptName
            print submitToQueue
            os.system(submitToQueue)
            time.sleep( 1.0 )

print "@@@@@ END JOB @@@@@@@@@@@@@@@"

