#!/bin/bash 
which grid-proxy-info
 source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
 echo $SCRAM_ARCH
 cmsenv
 source /cvmfs/cms.cern.ch/crab3/crab.sh
 which crab
 voms-proxy-init --voms cms --valid 168:00
