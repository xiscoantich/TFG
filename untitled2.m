%Script for Data.m
folder = fileparts(which(mfilename)); 
addpath(genpath(folder)); 

clearvars;
close all;
s.name = 'chirp';
s.type = 'TEMPORAL';
originalData = Data(s);

s.signal = originalData;
s.keep   = 0.5;
c = Compressor(s);
compressedData = c.computeCompressedSignal();

plot(originalData.signal);hold on;
plot(compressedData.signal);