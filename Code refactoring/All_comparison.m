%All comparison

in.typesignal = 'image'; % image or audio
in.filename = 'cat2.jpg'; %name of the file

data1 = Data(in);

figure();
set(gcf,'position',[0,0,680,3*440])
t1 = tiledlayout(3,3,'TileSpacing','compact','Padding','tight');

for i = [0.15 0.1 0.05]
c.keep = i;
c.method = 'threshold';

%Matlab
in.transtype = 'Fourier';
in.transmethod = 'dct_8by8'; 
t1 = Transformer(data1, in);
c_1 = Compressor(t1,c);
c_1.computeErr(data1);

%Packet
in.transtype = 'Wavelet';
in.transmethod = 'packet'; 
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
t2 = Transformer(data1, in);
c_2 = Compressor(t2,c);
c_2.computeErr(data1);

%PCA
in.transtype = 'PCA';
in.transmethod = ''; 
t3 = Transformer(data1, in);
c_3 = Compressor(t3,c);
c_3.computeErr(data1);



%% Plot

%colorbar

nexttile;
imshow(mat2gray(real(c_1.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['JPEG: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_1.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_2.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['Packet decom: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_2.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_3.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1

title(['PCA: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_3.err.mssim)],'FontSize',7);
end

%% 2 Error plot  
in.typesignal = 'image'; % image or audio
in.filename = 'cat.jpg'; %name of the file
data1 = Data(in);

N=50;
y1 = linspace(0.5,0.05,N);
in.transtype = 'Fourier';
in.transmethod = 'dct_8by8'; 
t1 = Transformer(data1, in);
in.transtype = 'Wavelet';
in.transmethod = 'packet'; 
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
t2 = Transformer(data1, in);
in.transtype = 'PCA';
in.transmethod = ''; 
t3 = Transformer(data1, in);

for i=1:N
c.keep = y1(i);
c.method = 'threshold';

%Matlab
c_1 = Compressor(t1,c);
c_1.computeErr(data1);

%Packet
c_2 = Compressor(t2,c);
c_2.computeErr(data1);

%PCA
c_3 = Compressor(t3,c);
c_3.computeErr(data1);

err1.ft(i) = c_1.err.mssim;
err1.wt(i) = c_2.err.mssim;
err1.pca(i) = c_3.err.mssim;
end

in.typesignal = 'image'; % image or audio
in.filename = 'cat2.jpg'; %name of the file
data1 = Data(in);

y1 = linspace(0.5,0.05,N);
in.transtype = 'Fourier';
in.transmethod = 'dct_8by8'; 
t1 = Transformer(data1, in);
in.transtype = 'Wavelet';
in.transmethod = 'packet'; 
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
t2 = Transformer(data1, in);
in.transtype = 'PCA';
in.transmethod = ''; 
t3 = Transformer(data1, in);


for i=1:N
c.keep = y1(i);
c.method = 'threshold';

%Matlab
c_1 = Compressor(t1,c);
c_1.computeErr(data1);

%Packet
c_2 = Compressor(t2,c);
c_2.computeErr(data1);

%PCA
c_3 = Compressor(t3,c);
c_3.computeErr(data1);

err2.ft(i) = c_1.err.mssim;
err2.wt(i) = c_2.err.mssim;
err2.pca(i) = c_3.err.mssim;
end


figure
t1 = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile
plot (y1*100,err1.ft,'DisplayName','JPEG');
hold on;
plot (y1*100,err1.wt,'DisplayName','Packet decom');
plot (y1*100,err1.pca,'DisplayName','PCA');
title('Image a)');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off

nexttile
plot (y1*100,err2.ft,'DisplayName','JPEG');
hold on;
plot (y1*100,err2.wt,'DisplayName','Packet decom');
plot (y1*100,err2.pca,'DisplayName','PCA');
title('Image b)');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off

%% 
%All comparison

in.typesignal = 'image'; % image or audio
in.filename = 'cat.jpg'; %name of the file

data1 = Data(in);

figure();
set(gcf,'position',[0,0,680,3*440])
t1 = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');

c1.keep = 0.02;
c1.method = 'threshold';

%Matlab
in.transtype = 'Fourier';
in.transmethod = 'dct_8by8'; 
t1 = Transformer(data1, in);
c_1 = Compressor(t1,c1);
c_1.computeErr(data1);

%Packet
c2.keep = 0.01;
c2.method = 'threshold';
in.transtype = 'Wavelet';
in.transmethod = 'packet'; 
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
t2 = Transformer(data1, in);
c_2 = Compressor(t2,c2);
c_2.computeErr(data1);

% %PCA
% in.transtype = 'PCA';
% in.transmethod = ''; 
% t3 = Transformer(data1, in);
% c_3 = Compressor(t3,c);
% c_3.computeErr(data1);



% Plot

%colorbar

nexttile;
imshow(mat2gray(real(c_1.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['JPEG: [',num2str(c1.keep*100),'%] MSSIM =', num2str(c_1.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_2.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['Packet decom: [',num2str(c2.keep*100),'%] MSSIM =', num2str(c_2.err.mssim)],'FontSize',7);

% nexttile;
% imshow(mat2gray(real(c_3.rec)));
% colormap;
% %colorbar %Los valores de la leyenda estan escalados entre 0 i 1
%title(['PCA: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_3.err.mssim)],'FontSize',7);
