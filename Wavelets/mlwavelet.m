function [c,l] = mlwavelet(Y,n,wavelet)
%load the wavelet here
if ischar(wavelet)
    wvf = load_wavelet(wavelet);
else
    wvf = wavelet;
end

Y=double(Y);
[Dcols]=size(Y,1);
subc=Dcols/(2^n);%dimensions of the lowest subband
%check the number of decompositions
if (round(subc) ~= subc)
    %at the moment, only powers of two supported
    error('Illegal number of decompositions for a given matrix!');
end

% l = zeros(1,n+1);
% c = zeros(1,length(Y));
% 
% for i = 1:1:n
%     [a,d] = dwt_lifting1D(Y,wvf);
%     l(end+1-i)=length(d);
%     c(1:length(a)+length(d)) = [a,d];
%     Y=a;
% end
% l(1) = length(a);


% Initialization.
s = size(Y);
Y = Y(:).'; % row vector
c = [];
l = zeros(1,n+2,'like',real(Y([])));

l(end) = length(Y);
for k = 1:n
    [Y,d] = dwt_lifting1D(Y,wvf); % decomposition
    c     = [d c];            %#ok<AGROW> % store detail
    l(n+2-k) = length(d);     % store length
end

% Last approximation.
c = [Y c];
l(1) = length(Y);

if s(1)>1
    c = c.'; 
    l = l';
end

end
