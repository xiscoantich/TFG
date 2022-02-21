classdef FourierWaveletsComparator < handle
    
    properties (Access = public)
        
            Fs = 1e2;                      %Frecuencia de muestreo de la senyal
            time =2;                      %Duracion de la señal
            keep = [0.9 0.1];
            
        %Wavelet parameters
            dj = 0.25;                    %smaller number gives better resolution, default = 0.25;
            Jfac = 1;                     %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default.
            mother = 'MORLET';            %The choices are 'MORLET', 'PAUL', or 'DOG'
            param = 6;                    %wave number for morlet, see >> help wave_bases for more details
    end
    
    properties (Access = private, Constant)
        %Example 1 
            f1=1;                           %Frecuencia de la señal Hz
            
        %Example 2
            f2_1=2;                         %Frecuencia de la señal Hz del 1º seno de construccion
            f2_2=3;                         %Frecuencia de la señal Hz del 2º seno de construccion
            
        %Example 3
            
        %Example 4
            
        %Wavelets
            pad = 1;
  
    end
    
    properties (Access = private, Dependent)
        %Wavelets
            dt;
            so;
            N;
            j1; 
    end
    
    methods 
        function value = get.dt(obj)
           value = 1/obj.Fs;
        end
        
        function value = get.so(obj)
           value = obj.dt;
        end
        
        function value = get.N(obj)
           value = obj.time*obj.Fs;
        end
        
        function value = get.j1(obj)
           value = round(obj.Jfac*(log2(obj.N*obj.dt/obj.so))/obj.dj); %default: (log2(N*dt/so))/dj
        end
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
            xMean = mean(x);
            x = x - xMean;
            [xWT, period, scale, coi, dj, paramout, k] = contwt(x, dt, pad, dj, so, j, mother, param); % contwt(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
            xCut = obj.cut(xWT,keep);
            xRec = invcwt(xCut, mother, scale, paramout, k) + xMean; %reconstruct the original signal            
        end

        function xRec = computeFourierSignal(obj,x,keep)
            xMean = mean(x);
            x = x - xMean;
            xF  = FFTCT_matrix(x);
            xCut = obj.cut(xF,keep);
            xRec = IFFTCT(xCut) + xMean;
            xRec = real(xRec);
        end        

        function computeExample1(obj,f1,Fs,time,keep,dt,pad,dj,so, j1, mother, param)
            figure
            [x,t] = obj.signal1(f1, Fs, time);
            obj.plot_signal(t, x);hold on;
            
            %Fourier
            for i=1:1:length(keep)
                xRec = obj.computeFourierSignal(x,keep(i));      
                plot(t,xRec,'-+','DisplayName',strcat('FFT r= ',num2str(keep(i))));
                hold on;
            end
            
            %Wavelet
            for i=1:1:length(keep)
                xW = obj.computeWaveletSignal(x, dt, pad, dj, so, j1, mother, param,keep(i));
                plot(t,xW,'-o','DisplayName',strcat('WT r= ',num2str(keep(i))));
                hold on;
            end
            legend('show')
          end

        function computeExample2(obj,f1, f2, Fs, time,dt,keep, pad, dj, so, j1, mother, param)
            figure
            [x,~,~,t] = obj.signal2(f1, f2, Fs, time);
            obj.plot_signal(t,x);
            hold on;
            
            %Fourier
            for i=1:1:length(keep)
                xF = obj.computeFourierSignal(x,keep(i));      
                plot(t,xF,'-+','DisplayName',strcat('FFT r= ',num2str(keep(i))));
                hold on;
            end
            
            %Wavelet
            for i=1:1:length(keep)
                xW = obj.computeWaveletSignal(x, dt, pad, dj, so, j1, mother, param,keep(i));
                plot(t,xW,'-o','DisplayName',strcat('WT r= ',num2str(keep(i))));
                hold on;
            end
            legend('show')
        end
   
         function computeExample3(obj, Fs, time,dt,keep, pad, dj, so, j1, mother, param)
            figure
            [x,t] = obj.signal3(Fs,time);
            obj.plot_signal(t,x);
            hold on;
            
            %Fourier
            for i=1:1:length(keep)
                xF = obj.computeFourierSignal(x,keep(i));      
                plot(t,xF,'-+','DisplayName',strcat('FFT r= ',num2str(keep(i))));
                hold on;
            end
            
            %Wavelet
            for i=1:1:length(keep)
                xW = obj.computeWaveletSignal(x, dt, pad, dj, so, j1, mother, param,keep(i));
                plot(t,xW,'-o','DisplayName',strcat('WT r= ',num2str(keep(i))));
                hold on;
            end
            legend('show')
         end


        function computeExample4(obj, Fs, time,dt,keep, pad, dj, so, j1, mother, param)
            figure
            [x,t] = obj.signal4(Fs,time);
            obj.plot_signal(t,x);
            hold on;
            
            %Fourier
            for i=1:1:length(keep)
                xF = obj.computeFourierSignal(x,keep(i));      
                plot(t,xF,'-+','DisplayName',strcat('FFT r= ',num2str(keep(i))));
                hold on;
            end
            
            %Wavelet
            for i=1:1:length(keep)
                xW = obj.computeWaveletSignal(x, dt, pad, dj, so, j1, mother, param,keep(i));
                plot(t,xW,'-o','DisplayName',strcat('WT r= ',num2str(keep(i))));
                hold on;
            end
            legend('show')
        end

        function compare(obj)    
            folder = fileparts(which(mfilename)); 
            % Add that folder plus all subfolders to the path.
            addpath(genpath(folder));   
            
            obj.computeExample1(obj.f1,obj.Fs,obj.time,obj.keep,obj.dt,obj.pad,obj.dj,obj.so,obj.j1,obj.mother,obj.param);

            obj.computeExample2(obj.f2_1, obj.f2_2, obj.Fs, obj.time,obj.dt,obj.keep, obj.pad, obj.dj, obj.so, obj.j1, obj.mother, obj.param);
            
            obj.computeExample3(obj.Fs, obj.time,obj.dt,obj.keep, obj.pad, obj.dj, obj.so, obj.j1, obj.mother, obj.param);

            obj.computeExample4(obj.Fs, obj.time,obj.dt,obj.keep, obj.pad, obj.dj, obj.so, obj.j1, obj.mother, obj.param)
            
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
            x2_1 = (sin(2*pi*f1*t2));   %1º Seno
            x2_2 = (sin(2*pi*f2*t2));   %2º Seno
            x2 = x2_1 + x2_2;           %Suma de senos
        end
        
        function [y, t] = signal3(Fs, time)     
            t = 0:1/Fs:time-(1/Fs);
            y = zeros(1,length(t));
            for i=1:length(t)
                y(i) =  piecewise_signal3(t(i),time);
            end
            y = y - mean(y);
            
            function y = piecewise_signal3(x,time)   %Funcion a trozos
                if x <= 0
                    y = 0;
                elseif 0 < x && x <= 1*time/4
                    y = (1*time/4)*6*x;
                elseif 1*time/4 < x && x <= 2*time/4 
                    y = 10*(2*time/4)*(x-(1*time/4))^2 -4*(2*time/4)*(x-1*time/4) + 0;
                elseif 2*time/4 < x && x <= 3*time/4
                    y = 0.5;
                else
                    y = 0;
                end 
            end
        end
        
        function [x,t] = signal4(Fs, time) %(Frequencia me muestreo, tiempo de la senyal)
            t = 0:1/Fs:time-(1/Fs);
            t1 = t(1:end/2);            % Primera mitad del vector
            t2 = t(end/2+1:end);        % Segunda mitad del vector
            x=[zeros(1,length(t1)) ones(1,length(t2))];
        end
        
        function plot_signal (t,x)
            plot(t,x,'DisplayName','Original');
            xlabel ('time [s]');
            ylabel ('g(t)');
            ylim([min(x)-(abs(max(x))*0.1) 1.1*max(x)])
            xlim([0 t(end)])
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