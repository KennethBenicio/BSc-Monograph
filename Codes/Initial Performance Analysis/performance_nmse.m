clc;
clear all;
close all;
pkg load communications;

runs = 10;
P = [2];
N = [4];
snr_db = [0 5 10 15 20 25 30];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 15);
graf_op = ['-kx'; '-bs';'-r+';'-mo';'-y*';];

nmse = zeros(length(snr_db),length(N));
for nn = 1:length(N)
    N(nn)
    for ii = 1:length(snr_db)
        for jj = 1:runs
            [nmse_test] = teste(snr_db(ii),P,N(nn));
            nmse(ii,nn) = nmse(ii,nn) + nmse_test;
        endfor
    endfor
    txt = ['N = ', num2str(N(nn))];
    semilogy(snr_db,nmse(:,nn)/runs,graf_op(nn,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt); 
endfor

title('IRS Performance: NMSE x SNR')
xlabel('SNR(dB)')
ylabel('NMSE')
grid on;
legend("location", "northeast");
hold off;