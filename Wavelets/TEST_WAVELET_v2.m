clear all;
t=0:0.01:10;
Y = sin(t);
% [wave, period, scale, coi, dj,paramout, k]= contwt(ref,1);
Nscales = 100; %contwt actually computes Nscales + 1 number of scales.
mothercwt = 'MORLET';
motherinv = 'MORLET';
DT = 1; %time step

% enter defaults for other parameters as []
[wave, period, scale, coi, dj,paramout, k]= contwt(Y,1,[], [], [], Nscales, mothercwt );
Xrec = invcwt(wave, motherinv, scale, paramout, k);

close all;
figure;

plot(Y);
hold on;
plot(Xrec, 'r');