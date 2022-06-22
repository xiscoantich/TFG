%% Fourier comparison
in.typesignal = 'image'; % image or audio
in.filename = 'cat.jpg'; %name of the file
in.transtype = 'Fourier';
data1 = Data(in);

figure();
%set(gcf,'position',[0,0,680,3*440])
t1 = tiledlayout(3,3,'TileSpacing','compact','Padding','tight');
%nexttile;
% imshow(mat2gray(data1.signal));
% title(['Original'],'FontSize',10);

%Matlab
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);

%dct
in.transmethod = 'dct'; 
ft2 = Transformer(data1, in);

%dct_8by8
in.transmethod = 'dct_8by8'; 
ft3 = Transformer(data1, in);

for i = [0.15 0.1 0.05]
c.keep = i;
c.method = 'threshold';


c_1 = Compressor(ft1,c);
c_1.computeErr(data1);

c_2 = Compressor(ft2,c);
c_2.computeErr(data1);

c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

% Plot

%colorbar

nexttile;
imshow(mat2gray(real(c_1.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['FFT: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_1.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_2.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['DCT: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_2.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_3.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1

title(['JPEG: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_3.err.mssim)],'FontSize',7);
end

%% 2 Error plot  
in.typesignal = 'image'; % image or audio
in.filename = 'cat.jpg'; %name of the file
in.transtype = 'Fourier';
data1 = Data(in);
N=20;
y1 = linspace(0.5,0.05,N);

%Matlab
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);
%dct
in.transmethod = 'dct'; 
ft2 = Transformer(data1, in);
%dct_8by8
in.transmethod = 'dct_8by8'; 
ft3 = Transformer(data1, in);

for i=1:N
c.keep = y1(i);
c.method = 'threshold';


c_1 = Compressor(ft1,c);
c_1.computeErr(data1);


c_2 = Compressor(ft2,c);
c_2.computeErr(data1);


c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

err1.fft(i) = c_1.err.mssim;
err1.dct(i) = c_2.err.mssim;
err1.jpeg(i) = c_3.err.mssim;
end

in.typesignal = 'image'; % image or audio
in.filename = 'cat2.jpg'; %name of the file
in.transtype = 'Fourier';
data1 = Data(in);

%Matlab
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);
%dct
in.transmethod = 'dct'; 
ft2 = Transformer(data1, in);
%dct_8by8
in.transmethod = 'dct_8by8'; 
ft3 = Transformer(data1, in);

for i=1:N
c.keep = y1(i);
c.method = 'threshold';


c_1 = Compressor(ft1,c);
c_1.computeErr(data1);


c_2 = Compressor(ft2,c);
c_2.computeErr(data1);


c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

err2.fft(i) = c_1.err.mssim;
err2.dct(i) = c_2.err.mssim;
err2.jpeg(i) = c_3.err.mssim;
end
figure
t1 = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile
plot (y1*100,err1.fft,'DisplayName','FFT');
hold on;
plot (y1*100,err1.dct,'DisplayName','DCT');
plot (y1*100,err1.jpeg,'DisplayName','JPEG');
title('Image a)');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off

nexttile
plot (y1*100,err2.fft,'DisplayName','FFT');
hold on;
plot (y1*100,err2.dct,'DisplayName','DCT');
plot (y1*100,err2.jpeg,'DisplayName','JPEG');
title('Image b)');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off