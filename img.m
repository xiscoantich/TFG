classdef img
    %Class to save and compres an image
        %1st u have to create the object and then introduce the image on .B
        %The compression is set to 1%
    properties 
        B
    end
    properties(Constant)
        keep=0.01;
    end
    properties(Dependent)
         Bt
        Atlow
        Alow
    end
    
    methods
        
        function value = get.Bt(obj)
           value = obj.FFT2CT_img(double(obj.B));
        end
        
        function value = get.Atlow(obj)
           Btsort = sort(abs(obj.Bt(:)));
           thresh = Btsort(floor((1-obj.keep)*length(Btsort)));
           ind = abs(obj.Bt)>thresh;       %Find small index;
           value = obj.Bt.*ind;            %Theshold small indices
        end
        
        function value =get.Alow(obj)
            value=uint8(ifft2(obj.Atlow));   %Compressed image
        end
        
        function plotFFT(obj)
            Blog = log(abs(fftshift(obj.Bt))+1); 
            figure
            imagen=mat2gray(Blog);
            %% Grey scale
            %imshow(img,[]); %Show the image in grey scale 
            %% Color scale
            imagesc(imagen);
            colormap;
            %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
            %% Mesh plot
            obj.plot_mesh(imagen);
        end
        
        function plot_FFTcut(obj)
            Blog = log(abs(fftshift(obj.Atlow))+1); %Put FFT on a logscale and offset 1 (fftshift(X) rearranges a Fourier transform X by shifting the zero-frequency component to the center of the array)
            figure
            imagen=mat2gray(Blog);
            %% Grey scale
            %imshow(img,[]);
            %title(['',num2str(keep*100),'%'],'FontSize',10)
            %% Color scale
            imagesc(imagen);
            colormap;
            %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
            title(['',num2str(obj.keep*100),'%'],'FontSize',10)
            %% Mesh plot
            obj.plot_mesh(imagen);
            title(['',num2str(obj.keep*100),'%'],'FontSize',10)
        end
        
        function plot_uncompress(obj)
            figure
            imshow(obj.Alow)               %Plot reconstruction
            title(['',num2str(obj.keep*100),'%'],'FontSize',10)
        end
        
        function plot_mesh (~,imagen)
            figure
            surf(imagen(10:10:end,10:10:end))
        end
        
        function [X]= FFT2CT_img(obj,x)
            %FFT2 Cooley-Tukey Algorithm
            m=size(x,2); %columnas
            for i=1:m
                x1(:,i)=obj.FFTCT_img(x(:,i)); %FFT por columnas
            end
            n=size(x1,1); %filas
            for i=1:n
                X(i,:)=obj.FFTCT_img(x1(i,:)); %FFT por filas
            end
        end
        
        function [X]= FFTCT_img(obj,x)
            %DIVIDE & CONQUER
            %Obtenemos un vector columna
            %Cuidado no puedes hacer las traspuesta del vector obtenido porque
            %reordena como le da la gana :')
            N1=length(x);
            %Check if its power of 2
            bool_power2=false;
            while bool_power2==false
                [x,bool_power2]=obj.makepowerof2(x,bool_power2);
            end
            N=length(x);
            %Separate the x[N] into even and odd-indexed subsequences
            for r=0:(N/2-1)
                %Even
                n=2*r;
                xe(r+1)=x(n+1);
                %Odd
                n=2*r+1;
                xo(r+1)=x(n+1);
            end
            if N<=2
                X(1,1)=xe+xo;
                X(2,1)=xe-xo;
            else
                Xe=obj.FFTCT_img(xe);
                Xo=obj.FFTCT_img(xo);  
                for k=0:length(xe)-1
                    w=exp(-(2*pi*1i/N)*k); 
                    X(k+1,1)=Xe(k+1)+w*Xo(k+1);
                    X(k+1+(N/2),1)=Xe(k+1)-w*Xo(k+1);
                end
            end
        end
        
        function [x,bool_power2]=makepowerof2(~,x,bool_power2)
            N=length(x);
            if mod(log(N)/log(2),1)~=0
                x(N+1)=0;
            else
                bool_power2=true;
            end
        end

    end
end

