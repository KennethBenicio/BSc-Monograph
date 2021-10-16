clc;
clear all;
close all;

format long;
N = [32, 64, 128, 256, 512, 1024];

figure
hold on;

graf_op = ['-ko'; '-m*'; '-yx'; '-bs'; '-gd'; '-r+';];

Pt = 0.01;
for nn = 1:length(N)
    
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
      angle_zeta(ww) = -(angle(G(ww)) + angle(H(ww)));
    endfor
    
    Z = diag(exp(1j.*angle_zeta));
    Y = (G.' * Z * H);
    
    snr_db =[0, 5, 10, 15, 20, 25, 30];
    %Capacity = zeros(length(snr_db),1);
    Error_rate = zeros(length(snr_db),1);

    runs = 1000;
    for ii = 1:length(snr_db)
        
        x = randn(N(nn),1) + 1j*randn(N(nn),1);
        snr_linear = 10^(snr_db(ii)/10);   
        variance_noise = 1/snr_linear;
        
        for jj = 1:runs    
            noise = (sqrt(variance_noise)/sqrt(2))*randn(N(nn),1);
            y = (G.' * Z * H)*sqrt(Pt)*x + noise;
            y = y/((G.' * Z * H)*sqrt(Pt));
            error = abs(y .- x);
            error = sum(error > 0.001);
            Error_rate(ii) = Error_rate(ii) + error/N(nn);
            %Capacity(ii) = Capacity(ii) + log2(1 + abs(Y)^2/sum(abs(noise).^2));     
        endfor    
    endfor
    
    %Capacity = Capacity/runs;
    Error_rate = Error_rate/runs;
    Error_rate(Error_rate<0) = 0;
    txt = ['N = ', num2str(N(nn))];
    %plot(snr_db,Capacity,graf_op(nn,:),"markersize", 15,'DisplayName',txt);
    plot(snr_db,Error_rate,graf_op(nn,:),"markersize", 15,'DisplayName',txt);

endfor

title('Error Rate vs SNR(dB) for N = {32, 64, 128, 256, 512, 1024}')
xlabel('SNR(dB)')
ylabel('Error Rate')
grid on;
legend("location", "southwest");
hold off;