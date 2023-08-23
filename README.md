# MedicalHerbChef


## About MedicalHerbChef
MedicalHerbChef (https://mhc.bio.kyutech.ac.jp/) enables us to predict the therapeutic effects and mode of action of multiherbal medicines consisting of arbitrary combinations and proportions of multiple crude drugs. Prediction is performed by machine learning based on molecular data (e.g., constituent compounds and chemical structures), omics-scale data (e.g., compound–protein interactomes and metabolome of crude drugs), and clinical data (e.g., formulas and efficacy). MedicalHerbChef accepts any combination and proportions of crude drugs as the input.

![overview](https://github.com/ShimadaYuuki/MedicalHerbChef/assets/100404818/4d6c0331-717b-453f-8d22-22dd638ea0db)


## Requirements

python 3.0 or later

gmpy2==2.0.0 or later

numpy==1.20.0 or later 

pandas==1.2.1 or later


## Usage

We provide two analysis modes:

**Custom-made proportion of crude drugs**

```
python　MedicalHerbChef.py　Dataset　Language　crude_drug1 dosage1 crude_drug2 dosage2
```


**A preset proportion of crude drugs**


```
python　MedicalHerbChef.py　Dataset　Language Kampo_formula
```
