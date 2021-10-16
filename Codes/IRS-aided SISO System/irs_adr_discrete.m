clc;
clear all;
close all;

N = [64];

figure
hold on;
ax = gca();
set(ax, 'fontsize', 12);
graf_op1 = ['--bs'; '--rs';'--gs'];
graf_op2 = ['-ks';];

bits = [1, 2, 3];

x_pilots =  randi([0,3],1,1);
x_pilots = 1/sqrt(2) * qammod(x_pilots,4); % 4-QAM Pilot Signal

for bb = 1:length(bits)
  
  QL = 2**bits(bb);
  fases = -QL/2:QL/2;
  fases(1) = [];
  angle_space = size(fases); 
  for ii = 1:QL
    angle_space(ii) =  2*pi*fases(ii)/QL; 
  endfor

  runs = 100;
  variance_signal = 1;
  snr_db =[0, 5, 10, 15, 20, 25, 30];
  Capacity1 = zeros(length(snr_db),1);
  Capacity2 = zeros(length(snr_db),1);
  for nn = 1:length(N)
      for jj = 1:runs    
        %z = (b-a)*rand() + a
        angle_g = pi.*rand(N(nn),1) - pi/2;
        angle_h = pi.*rand(N(nn),1) - pi/2;
        
        g = 1;
        h = 1;
        g = randn(N(nn),1);
        h = randn(N(nn),1);
        
        G = g .* exp(1j.*angle_g);
        G = G./abs(G);
        H = h .* exp(1j.*angle_h);
        H = H./abs(H);
        
        c_direct_link = rand(1,1) + 1j*rand(1,1);
      
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
        Y1 = G.' * Z1 * H * x_pilots + c_direct_link;
        Y2 = G.' * Z2 * H * x_pilots + c_direct_link;
      
        for ii = 1:length(snr_db)
          snr_linear = 10^(snr_db(ii)/10);   
          variance_noise = variance_signal/snr_linear;
          noise = (sqrt(variance_noise)/sqrt(2))*randn(N(nn),1);
          Capacity1(ii) = Capacity1(ii) + log2(1 + abs(Y1)^2/sum(abs(noise).^2));
          Capacity2(ii) = Capacity2(ii) + log2(1 + abs(Y2)^2/sum(abs(noise).^2));      
        endfor
      endfor
      
      if bb == 1
        Capacity2 = Capacity2/runs;
        txt2 = ['Continuous'];
        plot(snr_db,Capacity2,graf_op2(nn,:),'markersize', 12,'linewidth', 3,'markersize', 12, 'DisplayName',txt2);  
        Capacity1 = Capacity1/runs;
        txt1 = [num2str(bits(bb)),' bit of quantization'];
        plot(snr_db,Capacity1,graf_op1(bb,:), 'markersize', 12,'linewidth', 3, 'markersize', 12, 'DisplayName',txt1);  
      
      else
        Capacity1 = Capacity1/runs;
        txt1 = [num2str(bits(bb)),' bits of quantization'];
        plot(snr_db,Capacity1,graf_op1(bb,:), 'markersize', 12,'linewidth', 3, 'markersize', 12, 'DisplayName',txt1); 
       endif
  endfor
endfor

hold off;

title(['Continuous vs Discrete phase-shift N = ',num2str(N(nn))])
hx = xlabel('SNR in dB','fontsize',12);
hy = ylabel('SE in bps/Hz','fontsize',12);
set (hx, "fontweight", "bold")
set (hy, "fontweight", "bold")
        
%axis([0 30 0 30]);

legend_copy = legend("location", "southeast");
set (legend_copy, "fontsize", 16);

grid on;

print -depsc 2_continuous_vs_discrete.eps
