function example
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));  
d.type = 'TEMPORAL';
d.name = 'sinus';
d.motherwave = 'MORLET';
d.type_ft = 'matrix';

a = Data(d);

e.keep = 0.2;
e.data = a;
c = Compressor(e);
aCompressed = c.computeCompressedSignal();

figure()
hold on
a.plotSignal();
aCompressed.plotSignal()

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

end


%%  Errores o cosas a mejorar
%1. Cuando se hace la stft tambien quiere hacer despues una inversa despues y no deberia
%2. guardar las reconstrucciones porque no se como hacer para que no se sobre escriban las reconstrucciones ft y la wave
