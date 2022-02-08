classdef FourierWaveletsComparator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = FourierWaveletsComparator()
            obj.init()
            obj.compare();
        end
        
        
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
            
        end

        function xRec = computeWaveletSignal(obj,x, dt, pad, dj, so, j, mother, param,keep)
            [x_wt, period, scale, coi, dj, paramout, k] = contwt(x, dt, pad, dj, so, j, mother, param); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
            xCut = obj.cut(x_wt,keep);
            xRec = invcwt(xCut, mother, scale, paramout, k); %reconstruct the original signal            
        end

        function xRec = computeFourierSignal(obj,x,keep)
            xF  = FFTCT_matrix(x);
            xCut = obj.cut(xF,keep);
            xRec = IFFTCT(xCut);
            xRec = real(xRec);
        end        

        function computeExample1(obj,f1,Fs1,time,keep,dt,pad,dj,so, j1, mother, param)
            % Example 1
            figure
            [x,t] = obj.signal1(f1, Fs1, time);
            % x1 = x1(1:509); %507 - 513
            % t1 = t1(1:509);
            obj.plot_signal_1(t, x);hold on;
            
            %Fourier
            xRec = obj.computeFourierSignal(x,keep);      
            plot(t,xRec,'-+');
            hold on;
            
            %Wavelet
            xW = obj.computeWaveletSignal(x, dt, pad, dj, so, j1, mother, param,keep);

            
            plot(t,xW,'-o')
          end

        function computeExample2(obj,f1, f2, Fs, time,dt,keep, pad_2, dj, so, j1, mother, param)
            figure
            [x,x_1,x_2,t] = obj.signal2(f1, f2, Fs, time);
            obj.plot_signal_2(x,x_1,x_2,t);
            
            %Fourier
            xF = obj.computeFourierSignal(x,keep);      
            plot(t,xF,'-+');
            hold on;
            
            %Wavelet
            xW = obj.computeWaveletSignal(x, dt, pad_2, dj, so, j1, mother, param,keep);
            
            plot(t,xW,'-o')            
        end
   
        function computeExample3(obj,time3,keep3)
            figure
            
            [x, t3] = obj.signal3(time3);
            Fs=length(x)/time3;
            
            %Datos para wavelet
            pad_3 = 1;
            dj = 0.25;                    %smaller number gives better resolution, default = 0.25;
            dt = 1/Fs;
            so = dt;                    %default
            Jfac = 1;                       %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default.
            N = time3*Fs;
            j1 =  round(Jfac*(log2(N*dt/so))/dj); %default: (log2(N*dt/so))/dj
            mother_3 = 'MORLET';
            param_3 = 6;                    %wave number for morlet, see >> help wave_bases for more details
            
            %FFT
            x3_fft = FFTCT_matrix(x);
            p=zeros(1,length(keep3));
            p(1)=plot(t3,x);
            grid on;
            hold on
            for i=1:1:length(keep3)
                x3_fft_cut = obj.cut(x3_fft,keep3(i));
                x3_fft_rec = IFFTCT(x3_fft_cut);
                p(i+1)=plot(t3, real(x3_fft_rec)+mean(x),'LineWidth',1); %Aqui estoy sumando mean(x3) para compensar el DC
                hold on
            end
            
            %Wavelet
            for i=1:1:length(keep3)
                [x3_wt, period_3, scale_3, coi_3, dj, paramout_3, k_3] = contwt(x, dt, pad_3, dj, so, j1, mother_3, param_3); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
                x3_wt_cut = obj.cut(x3_wt,keep3(i));
                x3_wt_rec = invcwt(x3_wt_cut, mother_3, scale_3, paramout_3, k_3); %reconstruct the original signal
                p(length(keep3)+i+1)=plot(t3, real(x3_wt_rec)+mean(x),'LineWidth',1); %Aqui estoy sumando mean(x3) para compensar el DC
                hold on
            end
            
            legend([p(1) p(2) p(3) p(4) p(5)],{'Original','FFT r = 10%','FFT r = 90%','Wavelet r = 10%','Wavelet r = 90%'})

        end

        function computeExample4(obj,x4,time4,keep,dt,pad,dj,so, j1, mother, param,t4)
           [ft_4,x4_dft]= obj.ft_dft(x4, time4);           %Dft sin trenzado
            
            obj.plot_signal_1(t4, x4);
            hold on
            
            x4_fft = FFTCT_matrix(x4);
            x4_fft_cut = obj.cut(x4_fft,keep);
            x4_fft_rec = IFFTCT(x4_fft_cut);
            plot(t4,real(x4_fft_rec)+mean(x4));hold on; %Aqui estoy sumando mean(x4) para compensar el DC
            
            
            %Wavelet
            [x4_wt, period_4, scale, coi_4, dj_4, paramout, k_4] = contwt(x4, dt, pad, dj, so, j1, mother, param); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
            x4_wt_cut = obj.cut(x4_wt,keep);
            x4_wt_rec = invcwt(x4_wt_cut, mother, scale, paramout, k_4); %reconstruct the original signal
            plot(t4,real(x4_wt_rec)+mean(x4)); %Aqui estoy sumando mean(x4) para compensar el DC
            
            legend('Original','FFT r=25%','Wavelets r=25%');
        end



        
        function compare(obj)
            
            folder = fileparts(which(mfilename)); 
            % Add that folder plus all subfolders to the path.
            addpath(genpath(folder));   
            
            
            
            %Example 1
            Fs1 = 1e2;                      %Frecuencia de muestreo de la senyal
            time1 =6;                      %Duracion de la señal
            f1=1;                           %Frecuencia de la onda de la señal
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
            
            obj.computeExample1(f1,Fs1,time1,keep1,dt_1,pad_1,dj_1, so_1, j1_1, mother_1, param_1);
            
            
            %Example 2 
            
                        
            Fs2 = 1e2;                      %Frecuencia de muestreo de la senyal
            time2 =10;                      %Duracion de la señal
            f2_1=2;                         %Frecuencia de la señal Hz
            f2_2=3;                         %Frecuencia de la señal Hz            

            pad_2 = 1;
            dj_2 = 0.25;                    %smaller number gives better resolution, default = 0.25;
            dt_2 = 1/Fs2;
            so_2 = dt_1;                    %default
            Jfac_2 = 1;                       %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default.
            N_2 = time2*Fs2;
            j1_2 =  round(Jfac_2*(log2(N_2*dt_2/so_2))/dj_2); %default: (log2(N*dt/so))/dj
            mother_2 = 'MORLET';
            param_2 = 6;                    %wave number for morlet, see >> help wave_bases for more details

            obj.computeExample2(f2_1, f2_2, Fs2, time2,dt_2,keep1, pad_2, dj_2, so_2, j1_2, mother_2, param_2);


            
            % Example 3
            time3=10;            
            keep3=[0.1 0.9];              
            obj.computeExample3(time3,keep3)
            
            
            % Exampel 4
            time4 = 10;
            Fs4 = 2e1;
            
            figure
            x4=[zeros(1,time4/2*Fs4) ones(1,time4/2*Fs4)];
            t4=linspace(0,time4,Fs4*time4);
            
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
            keep_4 = 0.25;
            
            
            %Para la poder hacer la reconstrucción
            obj.computeExample4(x4,time4,keep_4,dt_4,pad_4,dj_4,so_4, j1_4, mother_4, param_4,t4)
            
        end
        
    end
    
    methods (Access = private, Static)
        
        
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
    end
    
    
    
    
end