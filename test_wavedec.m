clear all, close all, clc

dt = 0.05;
t = 0:dt:1-dt;
w = 15;
signal(:,1) = sin(w*t);

[c,l] = wavedec(signal,2,'Haar');
approx = appcoef(c,l,'Haar');
[cd1,cd2] = detcoef(c,l,[1,2]);

rec = waverec(c,l,'Haar');