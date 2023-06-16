# Tesis_Yago

## Contents
* [General Information](#General-Information)
* [Bibliografía](#Bibliografía)
* [Circulars_Automation](#Circulars_Automation)
* [Limits_GRBs](#Limits_GRBs)


## General Information

* Proyect title: 'Multimessenger astronomy with very-high-energy neutrinos at the observatory Pierre Auger Observatory'
* Autor : Yago Lema Capeáns
* Affiliations : [Universidade de Santiago de Compostela](https://www.usc.gal/es) & [Instituto Galego de Física de Altas Enerxías](https://igfae.usc.es/igfae/es/)

## Bibliografía

General Bibliography for the thesis 

## Circulars_Automation

The Pierre Auger Observatory is specifically designed to detect ultra-high-energy (UHE) cosmic rays,
making it a premier facility for studying these high-energy particles. However, it also possesses exceptional
sensitivity to UHE neutrinos across a significant portion of the sky. Of special relevance and potential
scientific impact is the followup in UHE neutrinos of the Gravitational Wave (GW) events detected by the
network of LIGO-Virgo-KAGRA (LVK) interferometric detectors. The followup of ∼ 100 GW events
collected in the O1, O2 and O3 runs of the LVK network has been performed in the past, with no
neutrino candidates identified, and reply circulars to the GCN alerts issued by LVK have been written and
sent manually. The O4 run of LVK has started on May 24, 2023 and will last 18 months. A significant
increase in the number of alerts on GW events that the Pierre Auger observatory is receiving is already
taking place, receiving around one event per day. In order to respond to them effectively we have created
software to automate the writing of circulars of the followup in UHE neutrinos. 

* Listenning_kafka once launched, it waits for CGN alerts from LVC collaboration (and IceCube), 
it runs subprocesses that perform plots and computations for each event, and finally
it sends a Circular automatically

* Contour.py is intended to obtain the Contours of 90% confidence level from
a GW event.

* fov_Auger_GW_builder.py has the aim to obtain the plot of the contours at 90% Confidence
Level of the GW event, together with the Field of view of the Pierre Auger.

* CL_Coverage plots the percentage of the 90% CL region of a GW
in coincidence with the FOV of the Pierre Auger observatory 
during an entire day after the detection time.

* Limit_flux.py Calculates the fluence limit assuming no detection and a normal
operating status of the SD1500 of the Pierre Auger Observatory 

## Limits_GRBs

Since no neutrino candidates have been identified so far in Auger data, we report on an update to the
individual and stacking fluence limits from neutrinos from Gamma-ray Bursts (GRBs) with the Surface
Detector of the Pierre Auger observatory. The limits are obtained assuming the Waxman-Bahcall model of
prompt neutrino emission in the GRB fireball that, at EeV energies, predicts a neutrino energy spectrum
dNν/dEν ∝ E_ν^(−4). We have identified a sample of 685 GRBs between 1 Jan 2004 and 31 Dec 2021 in the
field-of-view of the Surface Detector array, accounting for both the Earth-Skimming and the Downward-going 
high-angle channels that are the most sensitive to neutrinos. As a result of extending the data period
and hence the number of GRBs in the sample, and of including the Earth-Skimming channel, we have
improved by several orders of magnitude the stacking limit with respect to that reported in, reducing the
difference between the limit and the theoretical prediction of Waxman-Bahcall to a factor ≃ 5.
