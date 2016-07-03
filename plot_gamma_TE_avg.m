

Gamma_fit = [.225 .203 .174 .175 .178 .140 .132]/0.3;
Gamma_box = [.185 .19 .166 .159 .147 .126 .125]/0.3;
Gamma_fov = [.244 .218 .234 .221 .188 .176 .160]/0.3;
% rcore_fit = [0.0271 0.0224 0.0225 0.022 0.0252 0.0220 0.0196]/0.3;
deploy = [1 2 3 4 5 6 7];
labels = {'Fully Closed' '25% Open' '50% Open' '75% Open' 'Fully Open, X/C = 0' 'Fully Open, X/C = 1' 'Fully Open, X/C = 2'};
% plotyy(deploy,Gamma_fit,deploy,rcore_fit)
plot(deploy,Gamma_fit,'k',deploy,Gamma_box,'b',deploy,Gamma_fov,'g')%,deploy,rcore_fit,'k--')
ylim([0 1]);
legend('Vortex Strength (Parameter Fit)','Vortex Strength (Integral Box)','Vortex Strength (Integral FOV)','Location','Best')
% legend('Vortex Strength (Parameter Fit)','Vortex Strength (Integral Box)','Vortex Strength (Integral FOV)','Core Radius (Parameter Fit)','Location','Best')
set(gca, 'XTick', 1:7, 'XTickLabel', labels);
xlabel('Stage of DRS Deployment')
% ylabel('$\displaystyle\frac{dW}{dN}*$','interpreter','latex')
ylabel('$\displaystyle\frac{\Gamma}{Vc}$','interpreter','latex')
% ylabel('Vortex Strength (m^2/s); Radius (m)')