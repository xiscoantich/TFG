%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by 
% Abdullah Al Muhit
% contact - almuhit@gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M S Islam, and M Othman, “VLSI Implementation of Discrete Wavelet Transform (DWT) for Image Compression”, in Proc. of The Second International Conference on Autonomous Robots and Agents, Palmerston North, New Zealand, pp. 391-395, 2004, ISBN 0-476-00994-4. [PDF]
% A A Muhit, M S Islam, and M Othman, “ Design Design and Analysis of Discrete Wavelet Transform (DWT) for Image Compression Using VHDL”, in Proc. of the International Conference on Parallel and Distributed Processing Techniques and Applications, PDPTA 2005, Volume 1. CSREA Press 2005, pp. 157-160, Las Vegas, Nevada, USA, 2005, ISBN 1-932415-58-0. [PDF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE  demo of 1D signal and 2D image upto 1 level ,no compression

close all; clear all;
map = gray(256);

n=1:256;
x=sin((pi/64)*n)+0.01*(n-100);

% get db filter length 6
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('bior5.5');
%figure(1);
%subplot(2,2,1); plot(Lo_D);
%subplot(2,2,2); plot(Hi_D);
%subplot(2,2,3); plot(Lo_R);
%subplot(2,2,4); plot(Hi_R);

% 1D dwt and idwt
g = dwt(x,Lo_R);
y = idwt(g, Lo_R);
figure(2);
subplot(3,1,1); plot(x);
subplot(3,1,2); plot(g);
subplot(3,1,3); plot(y);

% 2D dwt and idwt
imData=imread('barbara.png');
figure(3);
imshow(imData);

% first dwt
[N, M]=size(imData);
% dwt in row
dwt_row_image=zeros(N, M);
tmpData=zeros(1, M);
for i=1:N
    tmpData(1, 1:M)=imData(i, 1:M);
    tmpData(1, 1:M)=dwt(tmpData, Lo_R);
    dwt_row_image(i, 1:M)=tmpData(1, 1:M);
end

%showData=zeros(N, M);
%showData(1:N, 1:M/2)=wcodemat(dwt_row_image(1:N, 1:M/2),256,'mat'); %wcodemat is command to normalize the image range
%showData(1:N, M/2+1:M)=wcodemat(dwt_row_image(1:N, M/2+1:M),256,'mat');
figure(4);
imshow(dwt_row_image, map);


% dwt in column
tmpData=zeros(1, N);
dwt1_imData=zeros(N, M);
for i=1:M
    tmpData(1, 1:N)=dwt_row_image(1:N, i)';
    tmpData(1, 1:N)=dwt(tmpData, Lo_R);
    dwt1_imData(1:N, i)=tmpData(1, 1:N)';
end

%showData=zeros(N, M);
%showData(1:N/2, 1:M/2)=wcodemat(dwt1_imData(1:N/2, 1:M/2),256,'mat');
%showData(1:N/2, M/2+1:M)=wcodemat(dwt1_imData(1:N/2, M/2+1:M),256,'mat');
%showData(N/2+1:N, 1:M/2)=wcodemat(dwt1_imData(N/2+1:N, 1:M/2),256,'mat');
%showData(N/2+1:N, M/2+1:M)=wcodemat(dwt1_imData(N/2+1:N, M/2+1:M),256,'mat');

figure(5);
imshow(dwt1_imData, map);

% my code for idwt
%idwt in coloumn

tmpData1=zeros(1, N);
idwt1_imData=zeros(N, M);
for i=1:M
    tmpData1(1, 1:N)=dwt1_imData(1:N, i)';
    tmpData1(1, 1:N)=idwt(tmpData1, Lo_R);
    idwt1_imData(1:N, i)=tmpData1(1, 1:N)';
end

figure(6);
imshow(idwt1_imData,map);

% idwt in row

idwt_row_image=zeros(N, M);
tmpData1=zeros(1, M);
for i=1:N
    tmpData1(1, 1:M)=idwt1_imData(i, 1:M);
    tmpData1(1, 1:M)=idwt(tmpData1, Lo_R);
    idwt_row_image(i, 1:M)=tmpData1(1, 1:M);
end

figure(7);
imshow(idwt_row_image,map);