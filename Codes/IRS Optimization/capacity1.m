clc;
clear all;
close all;

N = [32, 64, 128, 256, 512, 1024];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 15);
graf_op = ['-ko'; '-m*'; '-yx'; '-bs'; '-gd'; '-r+';];

bits = 3;
angle_space = linspace(-180,180,2**bits + 1)*(pi/180);
angle_space((2**bits + 2)/2) = [];

runs = 100;
variance_signal = 1;
snr_db =[0, 5, 10, 15, 20, 25, 30];
Capacity = zeros(length(snr_db),1);
for nn = 1:length(N)
    for jj = 1:runs    
      %z = (b-a)*rand() + a
      angle_g = pi.*rand(N(nn),1) - pi/2;
      angle_h = pi.*rand(N(nn),1) - pi/2;

      g = randn(N(nn),1);
      h = randn(N(nn),1);

      G = g .* exp(1j.*angle_g);
      G = G./abs(G);
      H = h .* exp(1j.*angle_h);
      H = H./abs(H);
    
      angle_zeta = zeros(N(nn),1);
      for ww = 1:N(nn)
        k = -(angle(G(ww))/2 + angle(H(ww))/2);
        index = interp1(angle_space,1:length(angle_space),k,'nearest');
        angle_zeta(ww) = angle_space(index);  
      endfor
    
      Z = diag(exp(1j.*angle_zeta));
      %channel with optimal phases.
      Y = G.' * Z * H;
      
      for ii = 1:length(snr_db)
        snr_linear = 10^(snr_db(ii)/10);   
        variance_noise = variance_signal/snr_linear;
        noise = (sqrt(variance_noise)/sqrt(2))*randn(N(nn),1);
        Capacity(ii) = Capacity(ii) + log2(1 + abs(Y)^2/sum(abs(noise).^2));      
      endfor
    endfor
    Capacity = Capacity/runs;
    txt = ['N = ', num2str(N(nn))];
    plot(snr_db,Capacity,graf_op(nn,:),"markersize", 15,'linewidth', 3,'DisplayName',txt);   
endfor

title(['Capacity vs SNR(dB) Considerando ' num2str(bits) ' bits de Quantização'])
xlabel('SNR(dB)')
ylabel('Capacity(bps/Hz)')
grid on;
legend("location", "northwest");
hold off;