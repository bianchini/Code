#!/usr/bin/env python

import commands
import re
import os

def processAllBatch( jobName , extra, t, g, b , G, T ):

    f = open('job_'+jobName+extra+'.sh','w')
    f.write('#!/bin/bash\n\n')
    f.write('cd /shome/bianchi/tthbb13/CMSSW/src\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('cd /shome/bianchi/Studies/Code\n')
    f.write('\n')
    f.write('./bin/toy -s 0 -t '+str(t)+' -p '+' -g '+str(g)+' -b '+str(b)+' -v 0'+' -G '+str(G)+' -T '+str(T)+' -o Test_'+jobName+'_smear'+str(g)+'_btag'+str(b)+'_gen'+str(G)+'_test'+str(T)+extra+'.root'+'\n')
    f.close()

    os.system('chmod +x '+'job_'+jobName+extra+'.sh')
    submitToQueue = 'qsub -V -cwd -l h_vmem=1G -q all.q -N job_'+jobName+extra+' job_'+jobName+extra+'.sh'
    os.system(submitToQueue)
    print submitToQueue


###########################################################################     

# TopLep + 5 jets radiation
#for p in range(20):
#    processAllBatch( 'tL_ggggg',  '_p'+str(p), 25, 1,1,  10,5 )
#    processAllBatch( 'lm_gggggg', '_p'+str(p), 25, 1,1,  11,5 )   

# TopHad + TopHad + Higgs
#for p in range(25):
#    processAllBatch( 'tH_tH_hH', '_p'+str(p), 20, 1,2,  4,4 )
#    processAllBatch( 'tH_tH_hH', '_p'+str(p), 20, 1,2,  9,4 )
#    processAllBatch( 'tH_tH_hH', '_p'+str(p), 20, 1,1, 14,4 )  
    
# TopHad + TopLep + Higgs
#for p in range(40):
#    processAllBatch( 'tH_tL_Hh', '_p'+str(p), 50, 1,1,  2,2 )
#    processAllBatch( 'tH_tL_Hh', '_p'+str(p), 50, 1,1,  7,2 )
#    processAllBatch( 'tH_tL_Hh', '_p'+str(p), 50, 1,1, 12,2 )
    
# TopLep + TopLep + Higgs
#for p in range(10):
#    processAllBatch( 'tL_tL_Hh', '_p'+str(p), 50, 1,1,  3,3 )
#    processAllBatch( 'tL_tL_Hh', '_p'+str(p), 50, 1,1,  8,3 )
#    processAllBatch( 'tL_tL_Hh', '_p'+str(p), 50, 1,1, 13,3 )

# TopHad + TopLep + Higgs: all hypos 
for p in range(50):
    # l+dddddd
    #processAllBatch( 'tH_tL_Hh-Ha1' , '_p'+str(p), 70, 1,1,  11,6 )

    # l+ddddbb
    #processAllBatch( 'tH_tL_Hh-Ha2' , '_p'+str(p), 70, 1,1,  15,6 )

    # tH+tL+dd
    #processAllBatch( 'tH_tL_Hh-Hb1' , '_p'+str(p), 70, 1,1,  12,6 )

    # tH+tL+cc
    #processAllBatch( 'tH_tL_Hh-Hb2' , '_p'+str(p), 70, 1,1,  16,6 )

    # tH+tL+bb
    #processAllBatch( 'tH_tL_Hh-Hb3' , '_p'+str(p), 70, 1,1,   7,6 )

    # tH+tL+Hh
    #processAllBatch( 'tH_tL_Hh-Hs'  , '_p'+str(p), 70, 1,1,   2,6 )

    # tH+tL+Hh: FullReco
    #processAllBatch( 'tH_tL_Hh-FullReco'  , '_p'+str(p), 20, 1,1,   2,  7 )
    # tH+tL+Hh: Lost
    #processAllBatch( 'tH_tL_Hh-NotReco'  ,  '_p'+str(p), 20, 1,1,   17, 7 )


    # tH+tL+Hh: FullReco                                                                                                               
    processAllBatch( 'tH_tL_Hh-FullReco-6jets'  , '_p'+str(p), 40, 1,1,   2,  8 )  

    # tH+tL+Hh: NotReco                                                                                                               
    processAllBatch( 'tH_tL_Hh-NotReco-6jets'   , '_p'+str(p), 40, 1,1,  17,  8 )  

    # tH+tL+dd: FullReco                                                                                                               
    processAllBatch( 'tH_tL_dd-FullReco-6jets'  , '_p'+str(p), 40, 1,1,  12,  8 )

    # tH+tL+dd: NotReco
    processAllBatch( 'tH_tL_dd-NotReco-6jets'   , '_p'+str(p), 40, 1,1,  18,  8 )

    # tH+tL+bb: FullReco
    processAllBatch( 'tH_tL_bb-FullReco-6jets'  , '_p'+str(p), 40, 1,1,   7,  8 )

    # tH+tL+bb: NotReco
    processAllBatch( 'tH_tL_bb-NotReco-6jets'   , '_p'+str(p), 40, 1,1,  19,  8 )
