function Image_Compression
close all; clearvars;
%test_Image_Compression
B = read_image;
plot_mesh(B);
Bt=fft2(B); 
plotFFT(Bt);
keep=0.01; %Esto se podria sacar fuera
Atlow = cut(Bt, keep);
plot_FFTcut(Atlow,keep);
Alow = uncompress(Atlow);
plot_uncompress(Alow, keep);
end

function B = read_image
A = imread('Images/cat.jpg');
B=rgb2gray(A);
end

function plotFFT(Bt)
Blog = log(abs(fftshift(Bt))+1); 
figure
img=mat2gray(Blog);
%% Grey scale
%imshow(img,[]);
%% Color scale
imagesc(img);
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
%% Mesh plot
plot_mesh(img);
end

function Atlow = cut(Bt, keep)
Btsort = sort(abs(Bt(:)));
thresh = Btsort(floor((1-keep)*length(Btsort)));
ind = abs(Bt)>thresh;       %Find small index;
Atlow = Bt.*ind;            %Theshold small indices
end

function plot_FFTcut(Atlow,keep)
Blog = log(abs(fftshift(Atlow))+1); %Put FFT on a logscale and offset 1 (fftshift(X) rearranges a Fourier transform X by shifting the zero-frequency component to the center of the array)
figure
img=mat2gray(Blog);
%% Grey scale
%imshow(img,[]);
%title(['',num2str(keep*100),'%'],'FontSize',10)
%% Color scale
imagesc(img);
colormap;
%colorbar %Los valores de la leyenda estan escalados entre 0 i 1
title(['',num2str(keep*100),'%'],'FontSize',10)
%% Mesh plot
plot_mesh(img);
title(['',num2str(keep*100),'%'],'FontSize',10)
end

function Alow =uncompress(Atlow)
Alow=uint8(ifft2(Atlow));   %Compressed image
end
function plot_uncompress(Alow, keep)
figure
imshow(Alow)               %Plot reconstruction
title(['',num2str(keep*100),'%'],'FontSize',10)
end

function plot_mesh (Img)
figure
surf(Img(10:10:end,10:10:end))
end

%% FFT test
function test_Image_Compression
load('test.mat');
Btest = test_read_image(test);
Bttest = test_FTT2(Btest,test);
Atlowtest = test_cut(Bttest,0.01,test); %Keep = 0.01
test_uncompress(Atlowtest,test);
end
% save ('Btest.mat','B')
% save ('Bttest.mat','Bt')
% save ('Atlowtest.mat','Atlow')
% save ('Alowtest.mat','Alow')

function Btest= test_read_image(test)
Btest = read_image;
Bdata = test.Btest;
error = norm(double(Bdata.B)-double(Btest));
if error==0
  cprintf('green', 'test read image  pass \n');

else
  cprintf('err', 'test read image error \n');
end
end

function Bttest= test_FTT2(Btest,test)
Bttest = fft2(Btest);
Btdata = test.Bttest;
error = norm(double(Btdata.Bt)-double(Bttest));
if error==0
  cprintf('green', 'test FTT2  pass \n');

else
  cprintf('err', 'test FTT2 error \n');
end
end

function Atlowtest = test_cut(Bttest,keep,test)
Atlowtest = cut(Bttest, keep);
Atlowdata = test.Atlowtest;
error = norm(double(Atlowdata.Atlow)-double(Atlowtest));
if error==0
  cprintf('green', 'test cut  pass \n');

else
  cprintf('err', 'test cut error \n');
end
end

function test_uncompress(Atlowtest,test)
Alowtest = uncompress(Atlowtest);
Alowdata = test.Alowtest;
error = norm(double(Alowdata.Alow)-double(Alowtest));
if error==0
  cprintf('green', 'test cut  pass \n');

else
  cprintf('err', 'test cut error \n');
end
end

