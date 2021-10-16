clc;
clear all;
close all;
pkg load communications;

runs = 100;
P = [1 2 3 4];
N = [16];
snr_db = [0 10 20 30];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 15);
graf_op = ['-kx'; '-bs';'-r+';'-mo';'-y*';];

adr = zeros(length(snr_db),length(P));
for pp = 1:length(P)
    pp
    for ii = 1:length(snr_db)
        for jj = 1:runs
            [_,adr_test] = IRS(snr_db(ii),P(pp),N);
            adr(ii,pp) = adr(ii,pp) + adr_test;
        endfor
    endfor
    txt = ['P = ', num2str(P(pp))];
    semilogy(snr_db,adr(:,pp)/runs,graf_op(pp,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt); 
endfor

title('IRS Performance: ADR x SNR')
xlabel('SNR(dB)')
ylabel('ADR(dB)')
grid on;
legend("location", "southeast");
hold off;