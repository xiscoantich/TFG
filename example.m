function example
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));  
d.type = 'TEMPORAL';
d.name = 'train';
d.motherwave = 'MORLET';
a = Data(d);


e.keep = 0.5;
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
aCompressed.plotWave()



end