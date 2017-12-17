clc
clear

tol1 = importdata('tol1.mat');
tol3 = importdata('tol3.mat');
tol5 = importdata('tol5.mat');
tol7 = importdata('tol7.mat');
tol9 = importdata('tol9.mat');

figure()

semilogy(1:100, abs(tol1), 'LineWidth', 2)

hold on

semilogy(1:100, abs(tol3), 'LineWidth', 2)
semilogy(1:100, abs(tol5), 'LineWidth', 2)
semilogy(1:100, abs(tol7), 'LineWidth', 2)
semilogy(1:100, abs(tol9), 'LineWidth', 2)

xlabel('Eiguenfrequency number')
ylabel('Difference betweent inverse iteration and eig (abs)')

legend('1e-01','1e-03', '1e-05', '1e-07', '1e-09')