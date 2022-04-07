clear all;

%load('ref.txt');
Y = rand(100,100);
% [wave, period, scale, coi, dj,paramout, k]= contwt(ref,1);
Nscales = 0; %contwt actually computes Nscales + 1 number of scales.
mothercwt = 'PAUL';
motherinv = 'PAUL';
DT = 1; %time step

% enter defaults for other parameters as []
[wave, period, scale, coi, dj,paramout, k]= contwt2(Y,1,[],[],[],[],mothercwt,[]);

%Ahora la inversa de bebe deshacer tambien siguiendo 
Xrec = invcwt(wave, motherinv, scale, paramout, k);

close all;
figure;

plot(Y);
hold on;
plot(Xrec, 'r');