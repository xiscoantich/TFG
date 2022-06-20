%Audio compressor
%% Fourier comparison
in.typesignal = 'audio'; % image or audio
in.filename = 'sinus'; %name of the file
in.transtype = 'Wavelet';
data1 = Data(in);
in.winfo.motherwave = 'Haar';

%Matlab
%in.transmethod = 'cwt'; 
%ft1 = Transformer(data1, in);

%dct
%in.transmethod = 'convolution'; 
%ft2 = Transformer(data1, in);

%dct_8by8
in.transmethod = 'lifting'; 
ft3 = Transformer(data1, in);

i = 0.15;
c.keep = i;
c.method = 'threshold';

%c_1 = Compressor(ft1,c);
%c_1.computeErr(data1);

%c_2 = Compressor(ft2,c);
%c_2.computeErr(data1);

c_3 = Compressor(ft3,c);
%c_3.computeErr(data1);
