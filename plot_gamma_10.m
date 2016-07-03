

Gamma_fit = [.225 .181 .159 .171 .181 .182 .138];
Gamma_box = [.210 .147 .145 .148 .159 .136 .118];
Gamma_fov = [.271 .241 .215 .244 .456 .153 .134];
rcore_fit = [.026 .016 .031 .015 .023 .023 .019];
% deploy = [0 .25 .50 .75 1 1.25 1.5];
deploy = [1 2 3 4 5 6 7];
labels = {'Fully Closed' '25% Open' '50% Open' '75% Open' 'Fully Open, X/C = 0' 'Fully Open, X/C = 1' 'Fully Open, X/C = 2'};

plot(deploy,Gamma_fit,'k',deploy,Gamma_box,'b',deploy,Gamma_fov,'g',deploy,rcore_fit,'k--')
legend('Vortex Strength (Parameter Fit)','Vortex Strength (Integral Box)','Vortex Strength (Integral FOV)','Core Radius (Parameter Fit)','Location','Best')
set(gca, 'XTick', 1:7, 'XTickLabel', labels);
xlabel('Stage of DRS Deployment')
ylabel('Vortex Strength (m^2/s); Radius (m)')