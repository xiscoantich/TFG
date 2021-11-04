%3blue1brown 
clearvars 
clc 
close all
% Pure signal
Fs = 1e2;                       %Frecuencia de muestreo
time =8;                        %Duracion de la señal
t = 0:1/Fs:time-(1/Fs);
f=4;                            %Frecuencia de la señal Hz
g1 = '(sin(2*pi*f*t))';
x = eval(g1);
figure;
plot(t,x);
xlabel ('time [s]');
ylabel ('g(t)');
ylim([-1.5 1.5])
xlim([0 time])
grid on;

% Wraping
figure
%Plot wraping around the origin;
tiledlayout(4,4);
for wf = 3.3:0.1:4.8 %winding frequency (cycles per second)3.3:0.1:4.8
nexttile
g2 = 'exp(-i*2*pi*wf.*t).*(sin(2*pi*f*t))';
x2 = eval(g2);
plot(real(x2), imag(x2))
ylim([-1 1])
xlim([-1 1])
title(['f = ',num2str(wf),'cycles/s'])
grid on
end

% Grafico f vs g(f)
%freq de muestreo para evitar trenzado
Fm=1/(2*f);
%Fm=1/(4*f);
%Calc graf
wf_vec = 0:Fm:50;

j=length(wf_vec);
ft_1=zeros(1,length(wf_vec));
for a=1:1:j
    wf=wf_vec(a);
    x2 = eval(g2);
    ft_1(a)=sum(x2);
end
figure
hold on
plot(wf_vec,abs(ft_1),'b')
grid on;

%FDT red 
xdft = fft(x);
xdft = xdft(1:length(x)/2+1);
freq = 0:Fs/length(x):Fs/2;
plot(freq,abs(xdft),'r+');
xlim([0 wf_vec(end)])
xlabel('f [Hz]');

%Recontruccion
%Obtener y recortar la funcion
%Aqui obtenemos la funcion pero multiplicado por 400
g_somb=cut(ft_1,0.01);

%Calculo
%g_somb=ft_1;
n = length(g_somb);
wf_vec_rec=wf_vec;
%t=linspace(0,8,n);
j=length(t);
for a=1:1:j
    time_rec=t(a);
    x1=g_somb.*exp(1i*2*pi*time_rec.*wf_vec_rec);
    x1_rec(a)=sum(x1);
end
%x3_sum=ifft(g_somb);
figure
plot(t,real(x1_rec)/length(wf_vec))
grid on;




%% Chord
clearvars;

Fs2 = 1e3;                       %Frecuencia de muestreo
time =10;                        %Duracion de la señal
t = 0:1/Fs2:time-(1/Fs2);
f1=2;                            %Frecuencia de la señal Hz
f2=3;                            %Frecuencia de la señal Hz
g1 = '(sin(2*pi*f1*t))';
g2 = '(sin(2*pi*f2*t))';
x1 = eval(g1);
x2 = eval(g2);
xsum=x1+x2;
figure;
tiledlayout(3,1)
nexttile
plot(t,xsum);
ylabel ('Dyad');
grid on;
nexttile
plot(t,x1);
ylabel ('g1(t)');
grid on;
nexttile
plot(t,x2);
grid on;
xlabel ('time [s]');
ylabel ('g2(t)');

%% Wrapping chord
figure
tiledlayout(4,4);
for wf = 1.8:0.1:3.3 %winding frequency (cycles per second)3.3:0.1:4.8
nexttile
g5 = '(exp(-i*2*pi*wf.*t).*(xsum))';
x5 = eval(g5);
plot(real(x5), imag(x5))
ylim([-2 2])
xlim([-2 2])
title(['f = ',num2str(wf),'cycles/s'])
grid on
end

%% Grafico f vs g(f)
wf_vec = 0:0.005:10;
j=length(wf_vec);
for a=1:1:j
    wf=wf_vec(a);
    x_ft = eval(g5);
    sumatorio(a)=sum(x_ft);
end

figure
plot(wf_vec,abs(sumatorio))
grid on;
hold on
%FDT red 
xdft = fft(xsum);
xdft = xdft(1:length(xsum)/2+1);
freq = 0:Fs2/length(xsum):Fs2/2;
plot(freq,abs(xdft));
xlim([0 wf_vec(end)])
xlabel('f [Hz]');
hold off;

%% Square wave 
clearvars
% time = (0:0.01:1)';
% y = 0.2*sin(2*pi*50*time) + 0.5*sin(2*pi*120*time);
% yn = y + 0.3*randn(size(time));
% figure
% plot(time,yn)
% %plot(t(1:50),yn(1:50))
t1 = 0:0.004:2;
y1 = exp(t1)-1;

t2 = 2:0.01:3;
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

time = linspace(0,10,n);
figure
plot(time,y)
grid on;
xlabel ('time [s]');
ylabel ('g(t)');

% Wraping sqw
figure
%Plot wraping around the origin;
tiledlayout(4,4);
g2_sqw = 'exp(-i*2*pi*wf.*time).*y';
for wf = 0:0.3:4.5 %winding frequency (cycles per second)3.3:0.1:4.8
nexttile
x2_sqw = eval(g2_sqw);
plot(real(x2_sqw), imag(x2_sqw))
ylim([-1 1])
xlim([-1 1])
title(['f = ',num2str(wf),'cycles/s'])
grid on
end

% Grafico f vs g(f) sqw
wf_vec_sqw = 0:0.01:6;
j=length(wf_vec_sqw);

for a=1:1:j
    wf=wf_vec_sqw(a);
    x2_sqw = eval(g2_sqw);
    sumatorio_ran(a)=sum(x2_sqw);
end
Fs=length(y)/10;
figure
plot(wf_vec_sqw,abs(sumatorio_ran))
grid on;
hold on
%FDT red 
xdft = fft(y);
xdft = xdft(1:length(y)/2+1);
freq = 0:Fs/length(y):Fs/2;
plot(freq,abs(xdft));
xlim([0 wf_vec_sqw(end)])
xlabel('f [Hz]');
hold off;

%Recontruccion
%Obtener y recortar la funcion
%Aqui obtenemos la funcion pero multiplicado por 1.6043e+04
g_somb=cut(sumatorio_ran,0.5);

%Calculo
%g_somb=ft_1;
n = length(g_somb);
wf_vec_rec=wf_vec_sqw;
%t=linspace(0,8,n);
j=length(time);
for a=1:1:j
    time_rec=time(a);
    x1=g_somb.*exp(1i*2*pi*time_rec.*wf_vec_rec);
    x1_rec(a)=sum(x1);
end
%x3_sum=ifft(g_somb);
figure
plot(time,real(x1_rec)/1.6043e+04)
grid on;
hold on 
plot(time,y)


%Funciones
function Atlow = cut(Bt, keep)
Btsort = sort(abs(Bt(:)));
thresh = Btsort(floor((1-keep)*length(Btsort)));
ind = abs(Bt)>thresh;       %Find small index;
Atlow = Bt.*ind;            %Theshold small indices
end
