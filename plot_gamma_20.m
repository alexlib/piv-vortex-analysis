

Gamma_fit = [.173 .208 .111 .168 .120 .166];
Gamma_box = [.188 .178 .124  .172 .148 .151];
Gamma_fov = [.243 .216 .177  .18 .135 .149];
rcore_fit = [.021 .023 .02 .021 .012 .027];
% deploy = [0 .25 .50 .75 1 1.25 1.5];
deploy = [1 2 3  5 6 7];
labels = {'Fully Closed' '25% Open' '50% Open' '75% Open' 'Fully Open, X/C = 0' 'Fully Open, X/C = 1' 'Fully Open, X/C = 2'};

plot(deploy,Gamma_fit,'k',deploy,Gamma_box,'b',deploy,Gamma_fov,'g',deploy,rcore_fit,'k--')
legend('Vortex Strength (Parameter Fit)','Vortex Strength (Integral Box)','Vortex Strength (Integral FOV)','Core Radius (Parameter Fit)','Location','Best')
set(gca, 'XTick', 1:7, 'XTickLabel', labels);
xlabel('Stage of DRS Deployment')
ylabel('Vortex Strength (m^2/s); Radius (m)')