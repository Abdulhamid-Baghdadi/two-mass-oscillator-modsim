g = 9.80665;
m_1 = 0.1;
m_2 = 0.1;
m_3 = 0.1;
m_s = m_1/2 + m_3;
m_ges = m_1 + m_2 + m_3;
k = 10;
c = 10;
l_1 = 0.2;
k_phi = 40;
c_phi = 0.4;
gain_4 = 1/(l_1^2*(m_1/3 + m_3))

F1 = 2;
M1 = 30;    
omega_0 = 10;

time = tout(:,1);
y_simulink = simout(:,1);
y_simscape = simout(:,2);
phi_simulink = simout(:,3);
phi_simscape = simout(:,4);

plot(time, y_simulink,'.','MarkerSize',15)
hold on
plot(time, y_simscape,'-','linewidth',2)
legend('y-simulink','y-simscape','Location','southeast')
title({'y: simulink vs simscape','F_1(t) = 2*cos(10*t) N', 'M_2(t) = 30*cos(10*t) Nm'}) 
xlabel('t [s]') 
ylabel('y [m]') 
hold off
figure()

plot(time, phi_simulink,'.','MarkerSize',15)
hold on
plot(time, phi_simscape,'-','linewidth',2)
legend('phi-simulink','phi-simscape','Location','southeast')
title({'phi: simulink vs simscape','F_1(t) = 2*cos(10*t) N', 'M_2(t) = 30*cos(10*t) Nm'}) 
xlabel('t [s]') 
ylabel('phi [m]') 
hold off