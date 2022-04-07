% FFT benchmark.
% Eight custom FFT implementations are benchmarked against the MatLab's built-in
% fft.m function.

clc; clear; close all;

N = 2.^(6:20);
L = length(N);

time_lapse = zeros(10,L);  % Here we keep the time each FFT implementation spends
                                            % until task completion.

for k = 1:L

   a = randn(1,N(k)) + 1i*randn(1,N(k));

    tic

        A1 = fft(a);

    time_lapse(1,k) = toc; 

    tic

        A2 = fft_rec(a);

    time_lapse(2,k) = toc; 

    tic

        A3 = fft_it(a);

    time_lapse(3,k) = toc; 

    tic

        A4 = radix2fft(a);

    time_lapse(4,k) = toc; 

    tic

        A5 = radix4fft(a);

    time_lapse(5,k) = toc; 

    tic

        A6 = fft_it2_mex(a);

    time_lapse(6,k) = toc; 
    
    tic

        A7 = splitradixfft(a);

    time_lapse(7,k) = toc; 
    
     tic

        A8 = dif_fft(a);

    time_lapse(8,k) = toc; 
    
    tic
    
        A9 = ger_fft(a);
        
    time_lapse(9,k) = toc;    
    
    tic
    
        A10 = ger_fft_mex(a);
        
    time_lapse(10,k) = toc;    

end 
    
%% Plot the performance
figure('Name','Benchmarking of FFT Performance');
semilogy(log2(N),time_lapse(1,:),'*-');  hold on;
semilogy(log2(N),time_lapse(2,:),'c*-');
semilogy(log2(N),time_lapse(3,:),'k*-');
semilogy(log2(N),time_lapse(4,:),'r*-');
semilogy(log2(N),time_lapse(5,:),'g*-');
semilogy(log2(N),time_lapse(6,:),'m*-');
semilogy(log2(N),time_lapse(7,:),'y*-');
semilogy(log2(N),time_lapse(8,:),'Color',[.5  .4  .7],'LineStyle','-','Marker','*');
semilogy(log2(N),time_lapse(9,:),'Color',[.4  .5  .7],'LineStyle','-','Marker','*');
semilogy(log2(N),time_lapse(10,:),'Color',[.2  .2  .6],'LineStyle','-','Marker','*');
grid;
axis tight;
legend('Built-In FFT','Radix-2, DIT, Recurcive FFT','Radix-2, DIT, Iterative FFT', ...
             'Radix-2, DIT, FFT', 'Radix-4, DIT, FFT','Radix-2, DIT, Iterative, mex-coded FFT',...
             'Split-Radix, DIT, FFT','Radix-2, DIF, FFT','Radix-2, DIT, "German" FFT', ...
             'Radix-2, DIT, "German" mex-coded FFT', 'Location','NorthWest');
title('Execution Time for Various FFT Implementations');
xlabel('log_2(length(x))');
ylabel('Execution Time (sec)');

%% Plot the Numerical Accuracy if you like !
% figure;
% subplot(7,2,1);
% stem(real(A1)-real(A2),'.')
% title('\Ree\{A1\} - \Ree\{A2\}');
% axis tight;
% grid on;
% 
% subplot(7,2,2);
% stem(imag(A1)-imag(A2),'r.')
% title('\Imm\{A1\} - \Imm\{A2\}');
% axis tight;
% grid on;
% 
% subplot(7,2,3);
% stem(real(A1)-real(A3),'.')
% title('\Ree\{A1\} - \Ree\{A3\}');
% axis tight;
% grid on;
% 
% subplot(7,2,4);
% stem(imag(A1)-imag(A3),'r.')
% title('\Imm\{A1\} - \Imm\{A3\}');
% axis tight;
% grid on;
% 
% subplot(7,2,5);
% stem(real(A1)-real(A4),'.')
% title('\Ree\{A1\} - \Ree\{A4\}');
% axis tight;
% grid on;
% 
% subplot(7,2,6);
% stem(imag(A1)-imag(A4),'r.')
% title('\Imm\{A1\} - \Imm\{A4\}');
% axis tight;
% grid on;
% 
% subplot(7,2,7);
% stem(real(A1)-real(A5),'.')
% title('\Ree\{A1\} - \Ree\{A5\}');
% axis tight;
% grid on;
% 
% subplot(7,2,8);
% stem(imag(A1)-imag(A5),'r.')
% title('\Imm\{A1\} - \Imm\{A5\}');
% axis tight;
% grid on;
% 
% subplot(7,2,9);
% stem(real(A1)-real(A6),'.')
% title('\Ree\{A1\} - \Ree\{A6\}');
% axis tight;
% grid on;
% 
% subplot(7,2,10);
% stem(imag(A1)-imag(A6),'r.')
% title('\Imm\{A1\} - \Imm\{A6\}');
% axis tight;
% grid on;
% 
% subplot(7,2,11);
% stem(real(A1)-real(A7),'.')
% title('\Ree\{A1\} - \Ree\{A7\}');
% axis tight;
% grid on;
% 
% subplot(7,2,12);
% stem(imag(A1)-imag(A7),'r.')
% title('\Imm\{A1\} - \Imm\{A7\}');
% axis tight;
% grid on;
% 
% subplot(7,2,13);
% stem(real(A1)-real(A8),'.')
% title('\Ree\{A1\} - \Ree\{A8\}');
% axis tight;
% grid on;
% 
% subplot(7,2,14);
% stem(imag(A1)-imag(A8),'r.')
% title('\Imm\{A1\} - \Imm\{A8\}');
% axis tight;
% grid on;