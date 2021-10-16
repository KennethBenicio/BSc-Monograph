clc;
clear all;
close all;

N = [128, 512, 1024];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 15);
graf_op1 = ['--kx'; '--bs';'--r+';];
graf_op2 = ['-kx'; '-bs';'-r+';];

bits = 3;
QL = 2**bits;
fases = -QL/2:QL/2;
fases(1) = [];
angle_space = size(fases); 
for ii = 1:QL
  angle_space(ii) =  2*pi*fases(ii)/QL; 
endfor

runs = 25;
variance_signal = 1;
snr_db =[0, 5, 10, 15, 20, 25, 30];
Capacity1 = zeros(length(snr_db),1);
Capacity2 = zeros(length(snr_db),1);
for nn = 1:length(N)
    for jj = 1:runs    
      %z = (b-a)*rand() + a
      angle_g = pi.*rand(N(nn),1) - pi/2;
      angle_h = pi.*rand(N(nn),1) - pi/2;

      %g = randn(N(nn),1);
      %h = randn(N(nn),1);
      g = 1;
      h = 1;

      G = g .* exp(1j.*angle_g);
      G = G./abs(G);
      H = h .* exp(1j.*angle_h);
      H = H./abs(H);
      
      angle_zeta = zeros(N(nn),1);
      for ww = 1:N(nn)
        k = -(angle(G(ww)) + angle(H(ww)));
        if sign(k) == 1
          [~,index] = min(abs(angle_space .- k));
          angle_zeta1(ww) = angle_space(index); 
        else 
          [~,index] = min(abs(k .- angle_space));
          angle_zeta1(ww) = angle_space(index); 
        endif
        
        angle_zeta2(ww) =  -(angle(G(ww)) + angle(H(ww)));
      endfor
    
      Z1 = diag(exp(1j.*angle_zeta1));
      Z2 = diag(exp(1j.*angle_zeta2));
      %channel with optimal phases.
      Y1 = G.' * Z1 * H;
      Y2 = G.' * Z2 * H;
      
      for ii = 1:length(snr_db)
        snr_linear = 10^(snr_db(ii)/10);   
        variance_noise = variance_signal/snr_linear;
        noise = (sqrt(variance_noise)/sqrt(2))*randn(N(nn),1);
        Capacity1(ii) = Capacity1(ii) + log2(1 + abs(Y1)^2/sum(abs(noise).^2));
        Capacity2(ii) = Capacity2(ii) + log2(1 + abs(Y2)^2/sum(abs(noise).^2));      
      endfor
    endfor
    Capacity1 = Capacity1/runs;
    Capacity2 = Capacity2/runs;
    txt1 = ['N = ', num2str(N(nn)),' Q'];
    txt2 = ['N = ', num2str(N(nn))];
    plot(snr_db,Capacity1,graf_op1(nn,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt1); 
    plot(snr_db,Capacity2,graf_op2(nn,:),"markersize", 15,'linewidth', 3, 'DisplayName',txt2);   
endfor

title([num2str(bits) ' bits de Quantização vs Não Quantizado'])
xlabel('SNR(dB)')
ylabel('Capacity(bps/Hz)')
%axis([0 30 0 30]);
grid on;
legend("location", "northwest");
hold off;