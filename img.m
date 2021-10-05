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
           value = fft2(obj.B);
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
    end
end

