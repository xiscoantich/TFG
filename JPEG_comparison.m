%% JPEG cut methods

in.typesignal = 'image'; % image or audio
in.filename = 'cat2.jpg'; %name of the file
in.transtype = 'Fourier';
data1 = Data(in);

% figure();
% set(gcf,'position',[0,0,680,3*440])
% t1 = tiledlayout(3,3,'TileSpacing','compact','Padding','tight');

N=20;
y1 = linspace(0.5,0.05,N);
for i=1:N
c.keep = y1(i);
c.method = 'threshold';
%dct_8by8
in.transmethod = 'dct_8by8'; 
ft3 = Transformer(data1, in);
c_3 = Compressor(ft3,c);
c_3.computeErr(data1);

c.keep = y1(i);
c.method = 'mask';
%dct_8by8
in.transmethod = 'dct_8by8'; 
ft4 = Transformer(data1, in);
c_4 = Compressor(ft4,c);
c_4.computeErr(data1);

err.jpegt(i) = c_3.err.mssim;
err.jpegr(i) = c_4.err.mssim;
end

figure
plot (y1*100,err.jpegt,'DisplayName','Threshold');
hold on;
plot (y1*100,err.jpegr,'DisplayName','Region');
xlabel('Compression [%]')
ylabel('MSSIM')
legend
hold off