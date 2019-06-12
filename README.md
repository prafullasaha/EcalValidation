# EcalValidation
Before Start follow these bellow steps:
```
cmsrel CMSSW_10_2_14

cd CMSSW_10_2_14/src

cmsenv

mkdir Validation

cd CMSSW_10_2_14/src/Validation

git clone https://github.com/prafullasaha/EcalValidation EcalValidation

cd $CMSSW_BASE/src

scram b -j 3
```
