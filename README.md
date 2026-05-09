The code was mainly developed by Ouwen Peng and Mengtian Jin
Project Overview
This repository contains the models and data associated with our research on the pressure-dependent mass transport and reaction kinetics involved in the electrochemical reduction of carbon dioxide (CO2RR)
To elucidate the underlying mechanisms of high-pressure CO2RR, this project relies on two primary computational models: a Diffusion-Reaction Model and a Kinetic Model

Diffusion Model
The diffusion-reaction model built upon the Nernst-Planck equation to describe the transport of CO2 and ionic species to the cathode surface
The model is used to simulate the local CO2, CO32- and HCO3- concentration and pH on the electrode surface as functions of applied cathode potential and pressure. 
The bult concentration of these species should be calculated first
The result will be used in following Kinetic model

Kinetic Model
The model is designed to quantify the fractional surface coverages of competing intermediates as a function of CO2 pressure and cathode potential 
The model  predicts the partial current densities of the resulting products across pressure
