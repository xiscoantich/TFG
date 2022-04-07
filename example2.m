function example2
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));  
d.type = 'TEMPORAL';
d.name = 'cat.jpg';
d.motherwave = 'MORLET';
d.type_ft = 'matlab';

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