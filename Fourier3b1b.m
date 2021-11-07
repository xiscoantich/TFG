%3blue1brown 
clearvars 
clc 
close all

%% Datos
%Example 1
Fs1 = 1e2;                     %Frecuencia de muestreo de la senyal
time1 =8;                      %Duracion de la señal
f1=4;
wf1 = 3.3:0.1:4.8;             %Wraping f
Fm1=1/60;                      %freq de muestreo (para evitar trenzado Fm=1/(2*f)
fend1=10;                      %Frecuencia mas alta para ft

%Example 2
Fs2 = 1e3;                       %Frecuencia de muestreo de la senyal
time2 =10;                       %Duracion de la señal
f2_1=2;                          %Frecuencia de la señal Hz
f2_2=3;                          %Frecuencia de la señal Hz
wf2 = 1.8:0.1:3.3;
Fm2=1/60;                      %freq de muestreo (para evitar trenzado Fm=1/(2*f)
fend2=10;                      %Frecuencia mas alta para ft

%Example 3
time3=10;
wf3 = 0:0.3:4.5;
Fs3 = 1e3;                     %Frecuencia de muestreo de la senyal
Fm3=1/60;                      %freq de muestreo (para evitar trenzado Fm=1/(2*f)
fend3=10;                      %Frecuencia mas alta para ft

%% Example 1
%Frecuencia de la señal Hz
[x1,t1] = signal1(f1, Fs1, time1);
plot_signal_1 (x1, t1);

% Wraping
plot_wraping (x1, wf1, t1)

% Grafico f vs g(f)
ft_1 = ft_3b1b (x1, Fs1, Fm1, fend1, t1);

%% Example 2
[x2,x2_1,x2_2,t2] = signal2(f2_1, f2_2, Fs2, time2);
plot_signal_2 (x2,x2_1,x2_2,t2);

% Wraping
plot_wraping (x2, wf2, t2);

% Grafico f vs g(f)
ft_2 = ft_3b1b (x2, Fs2, Fm2, fend2, t2);

%% Example 3

[x3, t3] = signal3(time3);
plot_signal_1(x3,t3);

% Wraping
plot_wraping (x3, wf3, t3)

% Grafico f vs g(f)
Fs3=length(x3)/10;
ft_3 = ft_3b1b (x3, Fs3, Fm3, fend3, t3);

%Recontruccion
%Obtener y recortar la funcion
%Aqui obtenemos la funcion pero multiplicado por 1.6043e+04
keep=[0.2 0.5 0.95];
p=zeros(1,length(keep));
wf_vec = 0:Fm3:fend3;
figure
p(1)=plot(t3,x3);
grid on;
hold on
for i=1:1:length(keep)
[x3_rec] = reconstruction (ft_3 ,keep(i), wf_vec, t3, 9.6597e+03);
p(i+1)=plot(t3, real(x3_rec));
end
legend([p(1) p(2) p(3) p(4)],{'Original','keep = 20%','keep = 50%','keep = 95%'})
hold off



%Funciones
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

function [ft_1] = ft_3b1b (x, Fs, Fm, f_end, t)
wf = 0:Fm:f_end;
j=length(wf);
ft_1=zeros(1,length(wf));
for i=1:1:j
    x2 = exp(-1i*2*pi*wf(i).*t).*x;
    ft_1(i)=sum(x2);
end
figure
hold on
plot(wf,abs(ft_1),'b','LineWidth',2)
grid on;

%FDT red 
xdft = fft(x);
xdft = xdft(1:length(x)/2+1);
freq = 0:Fs/length(x):Fs/2;
plot(freq,abs(xdft),'r','LineWidth',1.5);
xlim([0 f_end])
xlabel('f [Hz]');
legend('My code ft','Matlab fft');
hold off;
end

function [x1_rec] = reconstruction (sumatorio_ran ,keep, wf_vec_rec, time, divide)
%Recontruccion
%Obtener y recortar la funcion
g_somb=cut(sumatorio_ran,keep);
j=length(time);
for a=1:1:j
    time_rec=time(a);
    x1=g_somb.*exp(1i*2*pi*time_rec.*wf_vec_rec);
    x1_rec(a)=sum(x1);
end
x1_rec=x1_rec/divide;
end
