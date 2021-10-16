clc;
clear all;
close all;
pkg load communications;

runs = 500;
P = [1];
N = [32 64 128 256];
snr_db = [0 10 20 30];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 25);
graf_op1 = ['-kx'; '-bs';'-r+';'-mo';'-y*';];
graf_op2 = [':kx'; ':bs';':r+';':mo';':y*';];

adr1 = zeros(length(snr_db),length(N));
adr2 = zeros(length(snr_db),length(N));
for nn = 1:length(N)
    nn
    for ii = 1:length(snr_db)
        for jj = 1:runs
            [_,adr_test1] = IRS(snr_db(ii),P,N(nn));
            [_,adr_test2] = IRS2(snr_db(ii),P,N(nn));
            adr1(ii,nn) = adr1(ii,nn) + adr_test1;
            adr2(ii,nn) = adr2(ii,nn) + adr_test2;
        endfor
    endfor
    txt1 = ['Com Precoder N = ', num2str(N(nn))];
    txt2 = ['Sem Precoder N = ', num2str(N(nn))];
    plot(snr_db,adr1(:,nn)/runs,graf_op1(nn,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt1); 
    plot(snr_db,adr2(:,nn)/runs,graf_op2(nn,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt2); 
endfor

save adr1.mat adr1
save adr2.mat adr2

title('IRS Performance: ADR x SNR')
xlabel('SNR(dB)')
ylabel('ADR')
grid on;
legend("location", "southeast");
hold off;