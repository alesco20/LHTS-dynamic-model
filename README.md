# Fast dynamic model for Latent Heat Thermal Storage (LHTS) units
This is a python model for simulating the discharging process of a shell-and-tube LHTS tank. This code was developed by Alessandro Colangelo for his PhD thesis at Politecnico di Torino (Italy).

This modelling approach was initially developed studying the dynamic behaviour of a 40-kWh shell-and-tube LHTS unit from a system perspective. This storage unit is made of 96 copper pipes with longitudinal fins immersed in a medium-temperature phase change material (PCM). The heat transfer fluid (HTF) of a typical hydronic heating system flows inside the pipes. The model was also validated against experimental data of this storage unit.

### In few words
The model needs only a-priori known geometrical and thermo-physical features of the shell-and-tube LHTS unit that the modeller intends to study.

The inputs of the model are:
1. the inlet heat transfer fluid (HTF) mass flow rate at each time step
2. the inlet HTF temperature at each time step

The outputs of the model are: 
1. the HTF outlet temperature at each time step
2. the thermal power released by the LHTS at each time step
3. the LHTS state of charge (SOC) at each time step


## Using the model
[...]


## Model description
A thorough description of the model is contained in Chapter 4 of the following PhD thesis: "Thermal Energy Storage Technologies - Fast modelling, realisation and experimental characterisation of innovative latent heat storage units for system integration".

A visual representation of the modelling approach is presented in the following figure.
![Visual_abstract](https://github.com/alesco20/LHTS-dynamic-model/assets/116569046/d6ae8f69-1cdf-47e4-9fe5-69858304dd57)


## License 
Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
