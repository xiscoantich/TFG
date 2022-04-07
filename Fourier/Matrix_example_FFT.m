%Code for a 4x4 FFT
n=3;
%m=magic(n);
m=A;
% m1=zeros(n,n);
% m2=zeros(n,n);
for i=1:n
    m1(:,i)=FFTCT(m(:,i));
end
for i=1:n
    m2(i,:)=FFTCT(m1(i,:));
end

