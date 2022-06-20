%script_refactoring
function sctipt_imgCompression
    in.typesignal = 'image'; % image or audio
    in.filename = 'cat.jpg'; %name of the file
    data1 = Data(in); 
 
    in.transtype = 'Fourier'; %Fourier, Wavelet or PCA
    in.transmethod = 'matlab'; %Fourier: 
    %1D: matlab, matrix, dft,
    %2D: matlab, matrix, dct_8by8, dct
    ftdata1 = Transformer(data1, in);
    
    c.keep = 0.1;
    c.method = 'threshold';
    c_ftdata1 = Compressor(ftdata1,c);
    c_ftdata1.computeErr(data1);
    
    
    in2.typesignal = 'image'; % image or audio
    in2.filename = 'cat.jpg'; %name of the file
    data2 = Data(in2); 
 
    in2.transtype = 'Wavelet'; %Fourier, Wavelet or PCA
    in2.transmethod = 'packet'; %Fourier:
     %1D: matlab, matrix, dft,
    %2D: matlab, matrix, dct_8by8, dct
    in2.winfo.motherwave = 'CDF_9x7';%''db5';%'Haar';%'Haar';%'CDF_9x7';
    %Depending on the type of method some wavelets works and others don't
    %Matlab methods  --> wfilters
    %Own methods --> Wavelet_Data
    
    wtdata1 = Transformer(data2, in2);
    
    c.keep = 0.2;
    c.method = 'threshold';
    c_wtdata1 = Compressor(wtdata1,c);
    c_wtdata1.computeErr(data2);
    
    in.typesignal = 'image'; % image or audio
    in.filename = 'cat.jpg'; %name of the file
    data1 = Data(in); 
 
    in.transtype = 'PCA'; %Fourier, Wavelet or PCA
    in.transmethod = ''; %Fourier: 
    %1D: matlab, matrix, dft,
    %2D: matlab, matrix, dct_8by8, dct
    pcadata1 = Transformer(data1, in);
    
    c.keep = 0.1;
    c.method = 'threshold';
    c_pcadata1 = Compressor(pcadata1,c);
    c_pcadata1.computeErr(data1);

end