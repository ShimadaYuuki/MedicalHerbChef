# MedicalHerbChef


## About MedicalHerbChef
MedicalHerbChef (https://mhc.bio.kyutech.ac.jp/) enables us to predict the therapeutic effects and mode of action of multiherbal medicines consisting of arbitrary combinations and proportions of multiple crude drugs. Prediction is performed by machine learning based on molecular data (e.g., constituent compounds and chemical structures), omics-scale data (e.g., compound–protein interactomes and metabolome of crude drugs), and clinical data (e.g., formulas and efficacy). MedicalHerbChef accepts any combination and proportions of crude drugs as the input.

![overview](https://github.com/ShimadaYuuki/MedicalHerbChef/assets/100404818/4e5e65d3-04f9-41cb-bbfa-cd3bbb929043)

## Requirements

python==3.8.8

gmpy2==2.1.2

numpy==1.21.2

pandas==1.2.1


## Usage

Custom-made formula

```
python　MedicalHerbChef.py　KEGG　English　 Custom　Ephedra_Herb 5.0 Glycyrrhiza 1.5
```

Existing Kampo formula

```
python　MedicalHerbChef.py　KEGG　English　 Kakkonto
```
