# MedicalHerbChef


## About MedicalHerbChef
MedicalHerbChef (https://mhc.bio.kyutech.ac.jp/) enables us to predict the therapeutic effects and mode of action of multiherbal medicines consisting of arbitrary combinations and proportions of multiple crude drugs. Prediction is performed by machine learning based on molecular data (e.g., constituent compounds and chemical structures), omics-scale data (e.g., compound–protein interactomes and metabolome of crude drugs), and clinical data (e.g., formulas and efficacy). MedicalHerbChef accepts any combination and proportions of crude drugs as the input.

![overview](https://github.com/ShimadaYuuki/MedicalHerbChef/assets/100404818/4d6c0331-717b-453f-8d22-22dd638ea0db)


## Requirements

python==3.8.8

gmpy2==2.1.2

numpy==1.21.2

pandas==1.2.1


## Usage

We provide two analysis modes:

**Custom-made proportion of crude drugs**

```
python　MedicalHerbChef.py　KEGG　English　 Custom　Ephedra_Herb 5.0 Glycyrrhiza 1.5
```


**A preset proportion of crude drugs**


```
python　MedicalHerbChef.py　KEGG　English　 Kakkonto
```
