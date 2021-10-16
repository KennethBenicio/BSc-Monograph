clc;
clear all;
close all;
pkg load communications;

runs = 500;
P = [1];
N = [8 16 32];
snr_db = [0 10 20 30];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 15);
graf_op = ['-kx'; '-bs';'-r+';'-mo';'-y*';];

adr = zeros(length(snr_db),length(N));
for nn = 1:length(N)
    nn
    for ii = 1:length(snr_db)
        for jj = 1:runs
            [_,adr_test] = IRS(snr_db(ii),P,N(nn));
            adr(ii,nn) = adr(ii,nn) + adr_test;
        endfor
    endfor
    txt = ['N = ', num2str(N(nn))];
    plot(snr_db,adr(:,nn)/runs,graf_op(nn,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt); 
endfor

title('IRS Performance: ADR x SNR')
xlabel('SNR(dB)')
ylabel('ADR')
grid on;
legend("location", "southeast");
hold off;