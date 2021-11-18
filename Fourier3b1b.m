%3blue1brown new with dft
clearvars 
clc 
close all

%% Datos
%Example 1
Fs1 = 1e2;                     %Frecuencia de muestreo de la senyal
time1 =10;                      %Duracion de la señal
f1=4;
wf1 = 3.3:0.1:4.8;             %Wraping f
fend1=(1/time1)*(Fs1*time1/2); %Frecuencia mas alta para ft

%Example 2
Fs2 = 1e2;                       %Frecuencia de muestreo de la senyal
time2 =10;                       %Duracion de la señal
f2_1=2;                          %Frecuencia de la señal Hz
f2_2=3;                          %Frecuencia de la señal Hz
wf2 = 1.8:0.1:3.3;
fend2=(1/time2)*(Fs2*time2/2);   %Frecuencia mas alta para ft

%Example 3
time3=10;
wf3 = 0:0.3:4.5;
keep=[0.2 0.5 0.99];              %Reconstruccion

%Example 4
time4 = 10;
Fs4 = 1e2;

%% Example 1
%Frecuencia de la señal Hz
[x1,t1] = signal1(f1, Fs1, time1);
plot_signal_1 (x1, t1);
% Wraping
plot_wraping (x1, wf1, t1)

% Grafico f vs g(f)
[ft_1_trenz]= ft_3b1b (x1, Fs1, fend1, t1); %Para obtener trenzado
[ft_1,x1_dft]= ft_dft(x1, time1);           %Dft sin trenzado

%Reconstruction 
figure
[x1_rec] = reconstruction_dft (ft_1, 0.7, t1);

%% Example 2
[x2,x2_1,x2_2,t2] = signal2(f2_1, f2_2, Fs2, time2);
plot_signal_2 (x2,x2_1,x2_2,t2);

% Wraping
plot_wraping (x2, wf2, t2);

% Grafico f vs g(f)
[ft_2_trenz]= ft_3b1b (x2, Fs2, fend2, t2); %Para obtener trenzado
[ft_2]= ft_dft(x2, time2);                  %Dft sin trenzado

%Reconstruction 
figure
[x2_rec] = reconstruction_dft (ft_2, 0.7, t2);


%% Example 3

[x3, t3] = signal3(time3);
Fs3=length(x3)/time3;
fend3=(1/time3)*(Fs3*time3/2); %Frecuencia mas alta para ft
%plot_signal_1(x3,t3);

% Wraping
plot_wraping (x3, wf3, t3)

% Grafico f vs g(f)
[ft_3_trenz]= ft_3b1b (x3, Fs3, fend3, t3); %Para obtener trenzado
[ft_3]= ft_dft(x3, time3);                  %Dft sin trenzado   

%Reconstruction 

%Obtener y recortar la funcion
%Aqui obtenemos la funcion pero multiplicado por 1.6043e+04
p=zeros(1,length(keep));
figure
p(1)=plot(t3,x3);
grid on;
hold on
for i=1:1:length(keep)
[x3_rec] = reconstruction_dft (ft_3, keep(i), t3);
p(i+1)=plot(t3, real(x3_rec),'LineWidth',1.5);
hold on
end
legend([p(1) p(2) p(3) p(4)],{'Original','keep = 20%','keep = 50%','keep = 95%'})
hold off

%% Exampel 4
x4=[zeros(1,time4/2*Fs4) ones(1,time4/2*Fs4)];
t4=linspace(0,time4,Fs4*time4);
fend4=(1/time4)*(Fs4*time4/2); %Frecuencia mas alta para ft

plot_signal_1 (x4, t4);
% Wraping
% plot_wraping (x4, wf4, t4)

% Grafico f vs g(f)
[ft_4_trenz]= ft_3b1b (x4, Fs4, fend4, t4); %Para obtener trenzado
[ft_4,x4_dft]= ft_dft(x4, time4);           %Dft sin trenzado  

%Reconstruction 
figure
[x4_rec] = reconstruction_dft (ft_4, 0.1, t4);
hold on 
plot(t4,x4)

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

function plot_signal_1 (x, t)
figure;
plot(t,x);
xlabel ('time [s]');
ylabel ('g(t)');
ylim([min(x)-(abs(min(x))*0.1) 1.1*max(x)])
xlim([0 t(end)])
grid on;
end

function plot_signal_2 (x2,x2_1,x2_2,t2)
figure;
tiledlayout(3,1)
nexttile
plot(t2,x2);
ylabel ('Dyad');
grid on;
nexttile
plot(t2,x2_1);
ylabel ('g1(t)');
grid on;
nexttile
plot(t2,x2_2);
grid on;
xlabel ('time [s]');
ylabel ('g2(t)');
end

function plot_wraping (x, wf, t)
% Wraping around the origin;
figure
tiledlayout(2,8);
    for i=1:1:length(wf) %winding frequency (cycles per second)3.3:0.1:4.8
    nexttile
    x2 = exp(-1i*2*pi*wf(i).*t).*x;
    plot(real(x2), imag(x2))
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
figure
hold on
plot(wf,abs(ft_1),'b','LineWidth',1.5)
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
legend('My code ft','Matlab fft','FontSize', 12);
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
%Aqui lo graficamos todo
% xdft_cut = xdft(1:length(x));
% freq = 0:Fs/(length(x)):Fs-Fs/length(x); %Esto da mal
figure
plot(freq,abs(xdft_cut),'r+','LineWidth',1.5);
hold on;

%Plot DFT blue
x_rec_cut=x_rec(1:floor(length(x)/2+1));
plot(freq,abs(x_rec_cut),'bo','LineWidth',1.5)
grid on;
xlabel('f [Hz]');
legend('My code ft','Matlab fft','FontSize', 12);
hold off

%Plot error
error=(x_rec-xdft);
figure
plot(abs(error));
ylabel('error');
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

