This script combines data, screens those data for quality(/inclusion),
% and then performs frequency (and control analyses for the potential 
% effect of eye movements) to compare SSVEP strength after an injection 
% of saline to SSVEP strength after an injection of ketamine
%
% Companion code for:
%
% N-methyl d-aspartate receptor hypofunction reduces steady state visual
% evoked potentials (2023)
% Alexander Schielke & Bart Krekelberg
% Center for Molecular and Behavioral Neuroscience
% Rutgers University - Newark 

%1. combine lfp recordings and edf files (~15 minutes)
    combineFiles

%2. screen files for blinks and artifacts and determine whether files
%   should be included in further analyses (~5 minutes)
    screenFiles
    
%3. perform frequency analyses and analyze eye movements (~30 minutes)
    processData

%4. use linear mixed effects models to test hypothesis that ketamine
%   reduces strength of steady state visual evoked potentials (15 minutes)
    analyzeData

%5. visualize example data, data quality, frequency spectra and results
    createFigures