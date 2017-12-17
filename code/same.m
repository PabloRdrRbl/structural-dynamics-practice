clear
clc

iq = importdata('implicit_sol.mat');
eq = importdata('explicit_sol.mat');

T = 0.5 * (1 / 0.4651);
tf = 30 * T;
h1 = 1e-2;
h2 = 1e-4;

time1 = 0:h1:tf;
time2 = 0:h2:tf;

figure(1)

title('Solution for node 6')

subplot(1,3,1)

hold on

xlabel('Time [s]') 
ylabel('Displacement in the x-direction [m]') 

plot(time1, iq(37, :), 'r--', 'LineWidth', 3)
plot(time2, eq(37, :), 'b', 'LineWidth', 1)

ylim([-5.2e-3, 5.2e-3])
xlim([0, tf])

subplot(1,3,2)

hold on

xlabel('Time [s]') 
ylabel('Displacement in the y-direction [m]') 

plot(time1, iq(38, :), 'r--', 'LineWidth', 3)
plot(time2, eq(38, :), 'b', 'LineWidth', 1)

ylim([-5.2e-3, 5.2e-3])
xlim([0, tf])

subplot(1,3,3)

hold on

xlabel('Time [s]') 
ylabel('Displacement in the z-direction [m]') 

plot(time1, iq(39, :), 'r--', 'LineWidth', 3)
plot(time2, eq(39, :), 'b', 'LineWidth', 1)

legend('Implicit', 'Explicit')

ylim([-5.2e-3, 5.2e-3])
xlim([0, tf])