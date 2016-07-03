

Gamma_fit = [.257 .241 .188 .2 .225 .167 .145];
Gamma_box = [.248 .215 .189 .197 .165 .138 .138];
Gamma_fov = [.315 .268 .262 .205 .185 .188 .164];
rcore_fit = [0.034 0.03 0.027 0.029 0.028 0.024 0.022];
% deploy = [0 .25 .50 .75 1 1.25 1.5];
deploy = [1 2 3 4 5 6 7];
labels = {'Fully Closed' '25% Open' '50% Open' '75% Open' 'Fully Open, X/C = 0' 'Fully Open, X/C = 1' 'Fully Open, X/C = 2'};

plot(deploy,Gamma_fit,'k',deploy,Gamma_box,'b',deploy,Gamma_fov,'g',deploy,rcore_fit,'k--')
legend('Vortex Strength (Parameter Fit)','Vortex Strength (Integral Box)','Vortex Strength (Integral FOV)','Core Radius (Parameter Fit)','Location','Best')
set(gca, 'XTick', 1:7, 'XTickLabel', labels);
xlabel('Stage of DRS Deployment')
ylabel('Vortex Strength (m^2/s); Radius (m)')