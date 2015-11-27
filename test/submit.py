#!/usr/bin/env python
import commands
import time
import re
import os

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

splits = []
for k in range(10):
    splits.append([3*k+1,3*k+3])
print splits


count = 0
for part in splits:
    count += 1
    for rnd in [0]:
        for option in [6]:

            scriptName = 'job_'+str(rnd)+'_'+str(option)+'_p'+str(count)+'.sh'
            jobName    = 'job_'+str(rnd)+'_'+str(option)+'_p'+str(count)
            f = open(scriptName,'w')
            f.write('#!/bin/bash\n\n')
            f.write('cd /shome/bianchi/VHbb_3/CMSSW_7_4_15/src/MEIntegratorStandalone/test\n')
            f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
            f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
            f.write('eval `scramv1 runtime -sh`\n')
            f.write('\n')    
            f.write('../bin/btag_analyzer TTbar '+str(rnd)+' test '+str(part[0])+' '+str(part[1])+' '+str(option)+'\n')
            f.close()
            os.system('chmod +x '+scriptName)
            
            submitToQueue = 'qsub -V -cwd -q all.q -N '+jobName+' '+scriptName
            print submitToQueue
            os.system(submitToQueue)
            time.sleep( 0.5 )

print "@@@@@ END JOB @@@@@@@@@@@@@@@"

