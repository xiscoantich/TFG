function mse_error_plot
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
d.domain = 'TEMPORAL';
d.typesignal = 'IMAGE';%'AUDIO';
d.name = 'cat';%'sinus';
d.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
d.type.ft = 'matlab';%'dft';
d.type.wt = 'packet';%'dyadic_decomp';
d.level=5;
d.par=struct('N',d.level,'pdep',0,'wvf',d.motherwave,'dec','greedy');
d.ent_par=struct('ent','shannon','opt',0);

a = Data(d);
%y1 = logspace(0.15,0.000001,100) - 1;
%y1 = logspace(1,5,100)/10^5;
N=50;
y1 = linspace(0.2,0.001,N);
msr.fourier = zeros(1,N);
msr.wavelet = zeros(1,N);
msr.pca = zeros(1,N);
%plot(y1)
for i=1:N
e.keep = y1(i);
e.data = a;
c = Compressor(e);
aCompressed = c.computeCompressedSignal();
msr.fourier(i) = aCompressed.Fourier.error;
msr.wavelet(i) = aCompressed.Wavelet.error;
msr.pca(i) = aCompressed.PCA.error;
end
figure
plot (y1*100,msr.fourier,'DisplayName','Fourier');
hold on;
plot (y1*100,msr.wavelet,'DisplayName','Wavelet');
%plot (y1*100,msr.pca,'DisplayName','PCA');
xlabel('keep[%]')
ylabel('mse')
legend
hold off
end