clear all;
t=0:0.01:13.5-1;
Y1 = sin(t);
t2=0:0.01:(13.5-1)*2+0.01;

Y2 = sin(8*t);
Y = [Y1 Y2];
% [wave, period, scale, coi, dj,paramout, k]= contwt(ref,1);
Nscales = 200; %contwt actually computes Nscales + 1 number of scales.
mothercwt = 'PAUL';
motherinv = 'PAUL';
DT = 0.01; %time step

% enter defaults for other parameters as []
[wave, period, scale, coi, dj,paramout, k]= contwt(Y,1,[], [], [], Nscales, mothercwt );
Xrec = invcwt(wave, motherinv, scale, paramout, k);


img = flipud(mat2gray(real(wave(5:45, 1:end))));
[rows, columns] = size(img);
newWidth = round(0.1 * columns);
stretchedImage = imresize(img, [5*rows 2*newWidth]);

figure
t1 = tiledlayout(3,1,'TileSpacing','compact','Padding','tight');
nexttile 
plot (t2,Y)
xlim([0,25])
set(gca,'XTick',[], 'YTick', [])

nexttile([2 1])
imshow(stretchedImage);



close all;
figure;

plot(Y);
hold on;
plot(Xrec, 'r');

figure();
set(gcf,'position',[0,0,440,680/1.5])
t1 = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
nexttile
plot(Y);
xlim([0,2500])
xlabel('Time')
ylabel('')
set(gca,'XTick',[], 'YTick', [])
grid off

hold on;
nexttile([3,1]);
s=surf(real(wave));
s.EdgeColor = 'none';
ylim([20 40])
grid off
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
xlabel('b')
ylabel('a')
zlabel('W(a,b)')

