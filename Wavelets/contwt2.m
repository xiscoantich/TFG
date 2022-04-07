function [wave,period,scale,coi, dj, paramout, k]= contwt2(Y,dt,pad,dj,s0,J1,mother,param)
%[wave,period,scale,coi, dj, paramout, k] = contwt(Y,dt,pad,dj,s0,J1,mother,param);
%WAVE(N,b,c)    - N = length(Y) 
               %- b = # of scales 
               %- c = 1 for colums, 2 for rows
               
m=size(Y,2); %colums
for i=1:m
    [wave1(:,:,1),period,scale,coi, dj, paramout, k]=contwt(Y(:,i),dt,pad,dj,s0,J1,mother,param); %FFT por columnas
end
n=size(wave1,1); %rows
for i=1:n
    [wave(:,:,2),period,scale,coi, dj, paramout, k]=contwt(wave1(i,:),dt,pad,dj,s0,J1,mother,param); %FFT por filas
end
end