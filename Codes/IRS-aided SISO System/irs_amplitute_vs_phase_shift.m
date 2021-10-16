clc;
close all;
clear all;

phi       = pi/2; 
beta_min  = 0.2;
theta_irs = linspace(-pi,pi,10000);

figure

hold on;

alpha     = 1;
beta_n_theta_n_1 = (1 - beta_min).*((sin(theta_irs .- phi) .+ 1)/(2)).**alpha .+ beta_min;

txt = ['\alpha = ' num2str(alpha)];
plot(theta_irs, beta_n_theta_n_1, '--r', 'linewidth', 3,'displayname', txt);

alpha     = 2;
beta_n_theta_n_2 = (1 - beta_min).*((sin(theta_irs .- phi) .+ 1)/(2)).**alpha .+ beta_min;

txt = ['\alpha = ' num2str(alpha)];
plot(theta_irs, beta_n_theta_n_2, '--b', 'linewidth', 3,'displayname', txt);

alpha     = 0;
beta_n_theta_n_0 = (1 - beta_min).*((sin(theta_irs .- phi) .+ 1)/(2)).**alpha .+ beta_min;

txt = ['Ideal Amplitude'];

plot(theta_irs, beta_n_theta_n_0, '--k', 'linewidth', 3,'displayname', txt);

title(" Practical Relation Amplitude-Phase-Shift", "fontsize", 12);
hx = xlabel('Phase-Shift \theta_{n}','fontsize',12);
hy = ylabel('Amplitude \beta_{n}(\theta_{n})','fontsize',12);
set (hx, "fontweight", "bold")
set (hy, "fontweight", "bold")

axis([-pi pi 0 1])

legend_copy = legend("location", "southeast");
set (legend_copy, "fontsize", 12);
set(gca,'XTick',-pi:pi/2:pi) 
set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}, "fontsize", 12)

grid on;
hold off;

save beta_n_theta_n_1.mat beta_n_theta_n_1
save beta_n_theta_n_2.mat beta_n_theta_n_2
save beta_n_theta_n_0.mat beta_n_theta_n_0

%print -depsc 2_amplitude_vs_phase_shift.eps