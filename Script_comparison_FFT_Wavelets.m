%Code that compares FFT and wavelets
clc;
close all
clearvars
%This code uses the functions:
    %-FFTCT_matrix
    %-IFFTCT
    addpath 'C:\Users\Xisco Antich\Documents\UPC\TFG\Matlab Code\Fourier'
    %-contwt
    %-invcwt
    addpath 'C:\Users\Xisco Antich\Documents\UPC\TFG\Matlab Code\Wavelets'
    
 %% Datos
%Example 1
Fs1 = 1e2;                      %Frecuencia de muestreo de la senyal
time1 =6;                      %Duracion de la señal
f1=1;                           %Frecuencia de la onda de la señal
fend1=(1/time1)*(Fs1*time1/2);  %Frecuencia mas alta para ft
keep1 = 0.2;
%Datos para wavelet
pad_1 = 1;
dj_1 = 0.25;                    %smaller number gives better resolution, default = 0.25;
dt_1 = 1/Fs1;
so_1 = dt_1;                    %default
Jfac_1 = 1;                       %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default. 
N_1 = time1*Fs1;
j1_1 =  round(Jfac_1*(log2(N_1*dt_1/so_1))/dj_1); %default: (log2(N*dt/so))/dj
mother_1 = 'MORLET';
param_1 = 6;                    %wave number for morlet, see >> help wave_bases for more details

%Example 2
Fs2 = 1e2;                      %Frecuencia de muestreo de la senyal
time2 =10;                      %Duracion de la señal
f2_1=2;                         %Frecuencia de la señal Hz
f2_2=3;                         %Frecuencia de la señal Hz
wf2 = 1.8:0.2:3.3;
fend2=(1/time2)*(Fs2*time2/2);  %Frecuencia mas alta para ft

%Datos para wavelet
pad_2 = 1;
dj_2 = 0.25;                    %smaller number gives better resolution, default = 0.25;
dt_2 = 1/Fs2;
so_2 = dt_1;                    %default
Jfac_2 = 1;                       %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default. 
N_2 = time2*Fs2;
j1_2 =  round(Jfac_2*(log2(N_2*dt_2/so_2))/dj_2); %default: (log2(N*dt/so))/dj
mother_2 = 'MORLET';
param_2 = 6;                    %wave number for morlet, see >> help wave_bases for more details

%Example 3
time3=10;
wf3 = 0:0.3:4.5;
keep3=[0.1 0.9];                 %Reconstruccion


%Example 4
time4 = 10;
Fs4 = 2e1;
wf4 = 0:0.2:4.5;

%% Example 1
figure
[x1,t1] = signal1(f1, Fs1, time1);
% x1 = x1(1:509); %507 - 513
% t1 = t1(1:509);
plot_signal_1(t1, x1);hold on;

N = length(makepowerof2(x1));
%Fourier
x1_fft = FFTCT_matrix(x1);
x1_fft_cut = cut(x1_fft,keep1);
x1_fft_rec = IFFTCT(x1_fft_cut);
plot(t1,real(x1_fft_rec),'-+');hold on;

%Wavelet
[x1_wt, period_1, scale_1, coi_1, dj_1, paramout_1, k_1] = contwt(x1, dt_1, pad_1, dj_1, so_1, j1_1, mother_1, param_1); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
x1_wt_cut = cut(x1_wt,keep1);
x1_wt_rec = invcwt(x1_wt_cut, mother_1, scale_1, paramout_1, k_1); %reconstruct the original signal 

dE = sum(abs(x1 - x1_wt_rec))/length(x1); %compute point wise mean error
plot(t1,x1_wt_rec,'-o')


%% Example 2

figure
[x2,x2_1,x2_2,t2] = signal2(f2_1, f2_2, Fs2, time2);
plot_signal_2 (x2,x2_1,x2_2,t2);

N = length(makepowerof2(x2));hold on;
%Fourier
x2_fft = FFTCT_matrix(x2);
x2_fft_cut = cut(x2_fft,keep1);
x2_fft_rec = IFFTCT(x2_fft_cut);
plot(t2,real(x2_fft_rec),'-+');hold on;

%Wavelet
[x2_wt, period_2, scale_2, coi_2, dj_2, paramout_2, k_2] = contwt(x2, dt_2, pad_2, dj_2, so_2, j1_2, mother_2, param_2); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
x2_wt_cut = cut(x2_wt,keep1);
x2_wt_rec = invcwt(x2_wt_cut, mother_2, scale_2, paramout_2, k_2); %reconstruct the original signal 

dE = sum(abs(x2 - x2_wt_rec))/length(x2); %compute point wise mean error
plot(t2,x2_wt_rec,'-o')


%% Example 3
figure

[x3, t3] = signal3(time3);
Fs3=length(x3)/time3;
fend3=(1/time3)*(Fs3*time3/2); %Frecuencia mas alta para ft

%Datos para wavelet
pad_3 = 1;
dj_3 = 0.25;                    %smaller number gives better resolution, default = 0.25;
dt_3 = 1/Fs3;
so_3 = dt_3;                    %default
Jfac_3 = 1;                       %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default. 
N_3 = time3*Fs3;
j1_3 =  round(Jfac_3*(log2(N_3*dt_3/so_3))/dj_3); %default: (log2(N*dt/so))/dj
mother_3 = 'MORLET';
param_3 = 6;                    %wave number for morlet, see >> help wave_bases for more details

%FFT
x3_fft = FFTCT_matrix(x3);
p=zeros(1,length(keep3));
p(1)=plot(t3,x3);
grid on;
hold on
for i=1:1:length(keep3)
x3_fft_cut = cut(x3_fft,keep3(i));
x3_fft_rec = IFFTCT(x3_fft_cut);
p(i+1)=plot(t3, real(x3_fft_rec)+mean(x3),'LineWidth',1); %Aqui estoy sumando mean(x3) para compensar el DC
hold on
end

%Wavelet
for i=1:1:length(keep3)
[x3_wt, period_3, scale_3, coi_3, dj_3, paramout_3, k_3] = contwt(x3, dt_3, pad_3, dj_3, so_3, j1_3, mother_3, param_3); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
x3_wt_cut = cut(x3_wt,keep3(i));
x3_wt_rec = invcwt(x3_wt_cut, mother_3, scale_3, paramout_3, k_3); %reconstruct the original signal 
p(length(keep3)+i+1)=plot(t3, real(x3_wt_rec)+mean(x3),'LineWidth',1); %Aqui estoy sumando mean(x3) para compensar el DC
hold on
end

legend([p(1) p(2) p(3) p(4) p(5)],{'Original','FFT r = 10%','FFT r = 90%','Wavelet r = 10%','Wavelet r = 90%'})


%% Exampel 4
figure
x4=[zeros(1,time4/2*Fs4) ones(1,time4/2*Fs4)];
t4=linspace(0,time4,Fs4*time4);
fend4=(1/time4)*(Fs4*time4/2); %Frecuencia mas alta para ft

%Datos para wavelet
pad_4 = 1;
dj_4 = 0.25;                    %smaller number gives better resolution, default = 0.25;
dt_4 = 1/Fs4;
so_4 = dt_4;                    %default
Jfac_4 = 1;                       %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default. 
N_4 = time4*Fs4;
j1_4 =  round(Jfac_4*(log2(N_4*dt_4/so_4))/dj_4); %default: (log2(N*dt/so))/dj
mother_4 = 'MORLET';
param_4 = 6;                    %wave number for morlet, see >> help wave_bases for more details



%Para la poder hacer la reconstrucción
[ft_4,x4_dft]= ft_dft(x4, time4);           %Dft sin trenzado  

plot_signal_1 (t4, x4);
hold on

x4_fft = FFTCT_matrix(x4);
x4_fft_cut = cut(x4_fft,0.25);
x4_fft_rec = IFFTCT(x4_fft_cut);
plot(t4,real(x4_fft_rec)+mean(x4));hold on; %Aqui estoy sumando mean(x4) para compensar el DC


%Wavelet
[x4_wt, period_4, scale_4, coi_4, dj_4, paramout_4, k_4] = contwt(x4, dt_4, pad_4, dj_4, so_4, j1_4, mother_4, param_4); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
x4_wt_cut = cut(x4_wt,0.25);
x4_wt_rec = invcwt(x4_wt_cut, mother_4, scale_4, paramout_4, k_4); %reconstruct the original signal 
plot(t4,real(x4_wt_rec)+mean(x4)); %Aqui estoy sumando mean(x4) para compensar el DC

legend('Original','FFT r=25%','Wavelets r=25%');


%% Funciones
function Atlow = cut(Bt, keep)
Btsort = sort(abs(Bt(:)));
thresh = Btsort(floor((1-keep)*length(Btsort)));
ind = abs(Bt)>thresh;       %Find small index;
Atlow = Bt.*ind;            %Theshold small indices
end

function [g1,t1] = signal1(f, Fs, time) %(Frecuencia de la señal Hz, frequencia me muestreo, tiempo de la senyal)
t1 = 0:1/Fs:time-(1/Fs);            
g1 = (sin(2*pi*f*t1));
end

function [x2,x2_1,x2_2,t2] = signal2(f1, f2, Fs, time) %(Frecuencia de la señal Hz, frequencia me muestreo, tiempo de la senyal)
t2 = 0:1/Fs:time-(1/Fs);  
x2_1 = (sin(2*pi*f1*t2));
x2_2 = (sin(2*pi*f2*t2));
x2=x2_1+x2_2;
end

function [y, t] = signal3(time)
t1 = 0:0.004:2;
y1 = exp(t1)-1;

y2 = exp(2)-0.01-1:-0.01:3;

t3=0:0.01:3;
y3=5+(t3-t3);

t4=0:0.01:3.3;
a = 1.5;
b = -5;
c = 5;
y4=a*t4.^2+b*t4+c;

t6=0:0.01:3;
y5=2*sin(pi*t6)+4.9;

t6=-8:0.01:1.35;
y6=exp(t6)+1;
y6=flip(y6);

t7=0:0.01:5;
y7=2.5+(t7-t7);

y = [y1 y2 y3 y4 y5 y6 y7]/5-0.5;
n = length(y);

t = linspace(0,time,n);
end

function plot_signal_1 (t,x)
plot(t,x);
xlabel ('time [s]');
ylabel ('g(t)');
ylim([min(x)-(abs(max(x))*0.1) 1.1*max(x)])
xlim([0 t(end)])
grid on;
end

function plot_signal_2 (x2,x2_1,x2_2,t2)
ax1=nexttile(1,[2 2]);
plot(t2,x2);
ylabel ('Dyad');
grid on;
% xticklabels(ax1,{})
% ax2=nexttile(13,[3 2]);
% plot(t2,x2_1);
% ylabel ('g1(t)');
% grid on;
% xticklabels(ax2,{})
% nexttile(25,[3 2])
% plot(t2,x2_2);
% grid on;
xlabel ('time [s]');
% ylabel ('g2(t)');
end

function plot_wraping (x, wf, t)
% Wraping around the origin;
%tiledlayout(2,8);
    for i=1:1:8 %length(wf) %winding frequency (cycles per second)3.3:0.1:4.8
    nexttile
    x2 = exp(-1i*2*pi*wf(i).*t).*x;
    plot(real(x2), imag(x2));
    ylim([-1.1*max(x) 1.1*max(x)])
    xlim([-1.1*max(x) 1.1*max(x)])
    title(['f = ',num2str(wf(i)),'cycles/s'])
    grid on
    end
end

function [ft_1, xdft] = ft_3b1b (x, Fs, f_end, t)
n=10; %Para obtener trenzado 
% n=1; %Para no obtener trenzado
wf=linspace(0,f_end,n*length(t)/2+1);
j=length(wf);
ft_1=zeros(1,length(wf));
for i=1:1:j
    x2 = exp(-1i*2*pi*wf(i).*t).*x;
    ft_1(i)=sum(x2);
end

plot(wf,abs(ft_1),'b','LineWidth',1.5)
hold on
grid on;

%FDT red 
xdft = fft(x);
% Aqui solo estamos graficando la mitad de los puntos
xdft_cut = xdft(1:floor(length(x)/2+1));
%freq = 0:Fs/length(x):Fs/2; %Esto da mal 
freq = Fs*(0:(length(x)/2))/length(x);
%Aqui lo graficamos todo
% xdft_cut = xdft(1:length(x));
% freq = 0:Fs/(length(x)):Fs-Fs/length(x); %Esto da mal
plot(freq,abs(xdft_cut),'r','LineWidth',1);
xlabel('f [Hz]');
legend('My code ft','Matlab fft');
hold off;
end

function [x_rec,xdft] = ft_dft(x,time)
N=length(x);
Fs=N/time;
x_sum=zeros(1,N);
x_rec=zeros(1,N);
for k=0:1:N-1
    for n=0:1:N-1
    x_sum(n+1) = x(n+1)*exp(-(2*pi*1i/N)*k*(n));
    end
    x_rec(k+1) = sum(x_sum);
end

%FDT red 
xdft = fft(x);
% Aqui solo estamos graficando la mitad de los puntos
xdft_cut = xdft(1:floor(length(x)/2+1));
%freq = 0:Fs/length(x):Fs/2; %Esto da mal 
freq = Fs*(0:(length(x)/2))/length(x);
% %Aqui lo graficamos todo
% % xdft_cut = xdft(1:length(x));
% % freq = 0:Fs/(length(x)):Fs-Fs/length(x); %Esto da mal
% figure
% plot(freq,abs(xdft_cut),'r+','LineWidth',1.5);
% hold on;

% %Plot DFT blue
% x_rec_cut=x_rec(1:floor(length(x)/2+1));
% plot(freq,abs(x_rec_cut),'bo','LineWidth',1.5)
% grid on;
% xlabel('f [Hz]');
% legend('My code ft','Matlab fft','FontSize', 12);
% hold off

% %Plot error
% error=(x_rec-xdft);
% figure
% plot(abs(error));
% ylabel('error');
end

function [x1_rec] = reconstruction (sumatorio_ran ,keep, fend, time) %Esta funcion no sirve para nada
% Recontruccion
% Obtener y recortar la funcion
% wf_vec_rec = 0:Fm:fend;
wf_vec_rec=linspace(0,fend,length(time));
g_somb=cut(sumatorio_ran,keep);
j=length(time);
for a=1:1:j
    time_rec=time(a);
    x1=g_somb.*exp(1i*2*pi*time_rec.*wf_vec_rec);
    x1_rec(a)=sum(x1);
end
x1_rec=x1_rec/(length(sumatorio_ran)-1);
end

function [x_rec] = reconstruction_dft (ft, keep, t)

N=length(ft);
g_somb=cut(ft,keep);
x_sum=zeros(1,N);
x_rec=zeros(1,N);
for n=0:1:N-1
    for k=0:1:N-1
    x_sum(k+1) = g_somb(k+1)*exp((2*pi*1i/N)*k*(n));
    end
    x_rec(n+1) = 1/N*sum(x_sum);
end
plot(t,real(x_rec));
xlabel ('time [s]');
ylabel ('g(t)');
grid on;

end

function [y] = makepowerof2(x)
N = length(x);
y = x;
while mod(log(N)/log(2),1)~=0 
    y(N+1) = 0;
    N = N+1;
end
end
