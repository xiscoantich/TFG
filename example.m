function example
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
d.type = 'TEMPORAL';
d.typesignal = 'AUDIO';
d.name = 'train';
d.motherwave = 'CDF_9x7';
d.type_ft = 'matlab';
d.type_wt = 'lifting';
a = Data(d);

e.keep = 0.15;
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
        aCompressed.plotFrequency()
        
        figure()
        hold on
        a.plotWave()
        title('Wavelet')
        
        figure ()
        aCompressed.plotWave()
        title('Compressed Wavelet')
        
        figure()
        a.plotSurfWave()
        
        figure()
        aCompressed.plotSurfWave()
        
        figure()
        a.plotPCAInfo()
        
    case 'IMAGE' 
        figure()
        imagesc(a.signal);
        colorbar
        
        figure()
        imagesc(real(aCompressed.rec_f));
        colormap;
        %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
        title(['Fourier:',num2str(e.keep*100),'%'],'FontSize',10);
        
        figure()
        surf(real(aCompressed.rec_f(10:10:end,10:10:end)));
        
        figure()
        imagesc(real(aCompressed.rec_w));
        colormap;
        %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
        title(['Wavelets:',num2str(e.keep*100),'%'],'FontSize',10);
end


%%  Errores o cosas a mejorar
%1. Cuando se hace la stft tambien quiere hacer despues una inversa despues y no deberia
%2. guardar las reconstrucciones porque no se como hacer para que no se sobre escriban las reconstrucciones ft y la wave
% Esto tengo que crear un nuevo struct (clase) para cada uno de los metodos de
% compresion donde se guarde su recuperacion. No guardarlos todos en la
% misma clase. (Pasar de guardar en rec_f a signal)