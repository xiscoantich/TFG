%% Wavelet comparison
in.typesignal = 'image'; % image or audio
in.filename = 'cat2.jpg'; %name of the file
in.transtype = 'Wavelet';
data1 = Data(in);

figure();
set(gcf,'position',[0,0,680,3*440])
t1 = tiledlayout(3,3,'TileSpacing','compact','Padding','tight');
%nexttile;
% imshow(mat2gray(data1.signal));
% title(['Original'],'FontSize',10);

for i = [0.25 0.2 0.15]
c.keep = i;
c.method = 'threshold';

%Matlab
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);
c_1 = Compressor(ft1,c);
c_1.computeErr(data1);

%dwt
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'dwt_matlab'; 
ft2 = Transformer(data1, in);
c_2 = Compressor(ft2,c);
c_2.computeErr(data1);

%packet
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'packet'; 
ft3 = Transformer(data1, in);
c_3 = Compressor(ft3,c);
c_3.computeErr(data1);



%% Plot

%colorbar

nexttile;
imshow(mat2gray(real(c_1.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['Matlab decom: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_1.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_2.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['DWT: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_2.err.mssim)],'FontSize',7);

nexttile;
imshow(mat2gray(real(c_3.rec)));
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['Packet decom: [',num2str(c.keep*100),'%] MSSIM =', num2str(c_3.err.mssim)],'FontSize',7);
end

%% Error

N=70;
y1 = linspace(0.5,0.01,N);

for i=1:N
c.keep = y1(i);
c.method = 'threshold';

%Matlab
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);
c_1 = Compressor(ft1,c);
c_1.computeErr(data1);

%dwt
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'dwt_matlab'; 
ft2 = Transformer(data1, in);
c_2 = Compressor(ft2,c);
c_2.computeErr(data1);

%packet
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'packet'; 
ft3 = Transformer(data1, in);
c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

err.matlab(i) = c_1.err.mssim;
err.dwt(i) = c_2.err.mssim;
err.packet(i) = c_3.err.mssim;
end

figure
plot (y1*100,err.matlab,'DisplayName','Matlab decom');
hold on;
plot (y1*100,err.dwt,'DisplayName','DWT');
plot (y1*100,err.packet,'DisplayName','Packet decom');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off

%% 2 Error plot  
in.typesignal = 'image'; % image or audio
in.filename = 'cat.jpg'; %name of the file
in.transtype = 'Wavelet';
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
data1 = Data(in);
N=50;
y1 = linspace(0.5,0.001,N);
%Matlab
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);
%dwt
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'dwt_matlab'; 
ft2 = Transformer(data1, in);
%packet
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'packet'; 
ft3 = Transformer(data1, in);

for i=1:N
c.keep = y1(i);
c.method = 'threshold';


c_1 = Compressor(ft1,c);
c_1.computeErr(data1);


%c_2 = Compressor(ft2,c);
%c_2.computeErr(data1);


c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

err1.matlab(i) = c_1.err.mssim;
%err1.dwt(i) = c_2.err.mssim;
err1.packet(i) = c_3.err.mssim;
end

in.typesignal = 'image'; % image or audio
in.filename = 'cat2.jpg'; %name of the file
in.transtype = 'Wavelet';
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
data1 = Data(in);
%Matlab
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'matlab'; 
ft1 = Transformer(data1, in);
%dwt
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'dwt_matlab'; 
ft2 = Transformer(data1, in);
%packet
in.winfo.motherwave = 'Haar';%''db5';%'Haar';%'Haar';%'CDF_9x7';
in.transmethod = 'packet'; 
ft3 = Transformer(data1, in);

for i=1:N
c.keep = y1(i);
c.method = 'threshold';


c_1 = Compressor(ft1,c);
c_1.computeErr(data1);


%c_2 = Compressor(ft2,c);
%c_2.computeErr(data1);


c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

err2.matlab(i) = c_1.err.mssim;
%err2.dwt(i) = c_2.err.mssim;
err2.packet(i) = c_3.err.mssim;
end


figure
t1 = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile
plot (y1*100,err1.matlab,'DisplayName','Matlab decom');
hold on;
%plot (y1*100,err1.dwt,'DisplayName','DWT');
plot (y1*100,err1.packet,'DisplayName','Packet decom');
title('Image a)');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off

nexttile
plot (y1*100,err2.matlab,'DisplayName','Matlab decom');
hold on;
%plot (y1*100,err2.dwt,'DisplayName','DWT');
plot (y1*100,err2.packet,'DisplayName','Packet decom');
title('Image b)');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off