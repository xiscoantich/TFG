function example
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
d.type = 'TEMPORAL';
d.typesignal = 'IMAGE';
d.name = 'cat';
d.motherwave = 'Haar';
d.type_ft = 'matlab';
d.type_wt = 'dwt';
a = Data(d);

e.keep = 0.2;
e.data = a;
c = Compressor(e);
aCompressed = c.computeCompressedSignal();


switch d.typesignal
    case 'AUDIO'
        figure()
        hold on
        a.plotSignal('Original');
        aCompressed.Fourier.plotSignal('Fourier')
        aCompressed.Wavelet.plotSignal('Wavelet')
        aCompressed.PCA.plotSignal('PCA')
        
        figure()
        hold on
        a.plotFrequency();
        aCompressed.Fourier.plotFrequency()
        
        figure()
        hold on
        a.plotWave()
        title('Wavelet')
        
        figure ()
        aCompressed.WaveletplotWave()
        title('Compressed Wavelet')
        
        figure()
        a.plotSurfWave()
        
        figure()
        aCompressed.Wavelet.plotSurfWave()
        
        figure()
        a.plotPCAInfo()
        
    case 'IMAGE' 
        
        figure()
        surf(real(aCompressed.Fourier.signal(10:10:end,10:10:end)));
        
        figure()
        t = tiledlayout(4,4,'TileSpacing','compact','Padding','tight');
        nexttile([2,2]);
        imshow(mat2gray(a.signal));
        title(['Original'],'FontSize',10);
        %colorbar
        
        nexttile([2,2]);
        imshow(mat2gray(real(aCompressed.Fourier.signal)));
        colormap;
        %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
        title(['Fourier:',num2str(e.keep*100),'%'],'FontSize',10);
        
        
        nexttile([2,2]);
        imshow(mat2gray(real(aCompressed.Wavelet.signal)));
        colormap;
        %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
        title(['Wavelets:',num2str(e.keep*100),'%'],'FontSize',10);
        
        nexttile([2,2]);
        imshow(mat2gray(real(aCompressed.PCA.signal)));
        colormap;
        %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
        title(['PCA:',num2str(e.keep*100),'%'],'FontSize',10);
        
        f=figure();
        f.Position(3:4) = [600 250];
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
        nexttile;
        Blog = log(abs(fftshift(real(aCompressed.Fourier.freq))+1)); %Put FFT on a logscale and offset 1 (fftshift(X) rearranges a Fourier transform X by shifting the zero-frequency component to the center of the array)
        img=mat2gray(Blog);
        imshow(img);
        title(['Fourier coefficients:',num2str(e.keep*100),'%'],'FontSize',10);
        
        nexttile;
        Coef_wt=[aCompressed.Wavelet.wave(:,:,1),aCompressed.Wavelet.wave(:,:,2);
                 aCompressed.Wavelet.wave(:,:,3),aCompressed.Wavelet.wave(:,:,4)];
        img=mat2gray(Coef_wt);
        imshow(img);
        title(['Wavelets coefficients:',num2str(e.keep*100),'%'],'FontSize',10);
               
end


%%  Errores o cosas a mejorar
%1. Cuando se hace la stft tambien quiere hacer despues una inversa despues y no deberia
%2. guardar las reconstrucciones porque no se como hacer para que no se sobre escriban las reconstrucciones ft y la wave
% Esto tengo que crear un nuevo struct (clase) para cada uno de los metodos de
% compresion donde se guarde su recuperacion. No guardarlos todos en la
% misma clase. (Pasar de guardar en rec_f a signal)