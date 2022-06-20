%Threshold for report

in.typesignal = 'image'; % image or audio
in.filename = 'test.jpg'; %name of the file
data1 = Data(in);
padd = 972-688;
data1.signal = imresize(data1.signal, [688, 688]);
data1.signal = padarray(data1.signal,[padd/2, padd/2],'symmetric','both');
%Fourier
in.transtype = 'Fourier';
in.transmethod = 'matlab'; 
t1 = Transformer(data1, in);
figure();
set(gcf,'position',[0,0,688*3/2.5,2*688/2.5])
t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
nexttile
imshow(mat2gray(real(data1.signal(padd/2+1:688+padd/2, padd/2+1:688+padd/2))))



for i = [10 30 60 160 460]
    npixel = round(pi*(i)^2);
    k = npixel/(size(data1.signal,1)*size(data1.signal,2));
    c.keep = k;
    c.method = 'threshold';
    %Fourier
    c_1 = Compressor(t1,c);
    c_1.computeErr(data1); 
    nexttile
    imshow(mat2gray(real(c_1.rec(padd/2+1:688+padd/2, padd/2+1:688+padd/2))));
end