# EcalValidation
* Before Start follow these steps:
```
cmsrel CMSSW_10_2_14
cd $CMSSW_BASE/src
cmsenv
mkdir Validation
cd $CMSSW_BASE/src/Validation
git clone https://github.com/prafullasaha/EcalValidation EcalValidation
cd $CMSSW_BASE/src
scram b -j 3
```
* To run PFRechits part, go to 
`$CMSSW_BASE/src/Validation/EcalValidation/test/ecalvalidationMCAOD_cfg.py` 
file and make sure `usePFRecHitFlag` is `True`
* For interactive run:
```
voms-proxy-init --rfc --voms cms --valid 168:00
cd $CMSSW_BASE/src/Validation/EcalValidation/test
cmsRun ecalvalidationMCAOD_cfg.py
```
