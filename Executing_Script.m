MU_Effective = Permeability_Function();
E_effective = Permittivity_Function();

figure
step_num = 0.00001:0.000005:1; 
plot(step_num,MU_Effective,'--r', step_num, E_effective,'b', 'Linewidth', 1.2)
yline(0,'-.k','Linewidth', 1.2)
legend('$\mu_{eff}$', '$\varepsilon_{eff}$','Interpreter','Latex')
xlabel('$\omega_0^2/ \omega_p^2$','Interpreter','Latex','FontSize', 15) 
xlim([0.117345 0.121522]) % this needs to change for a different geometry
ylim([-30 5]) % this needs to change for a different geometry