%Script for Data.m
clearvars;
close all;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));  
d.type = 'TEMPORAL';
d.name = 'train';

a = Data(d);

e.keep = 0.5;
e.signal = a.signal;
c = Compressor(e);
computeCompressedSignal(c);


