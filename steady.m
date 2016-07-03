figure;
% [AX,H1,H2] = plotyy(X_over_C,[avg_Gamma_fit_vect; avg_Gamma_box_vect; avg_Gamma_FOV_vect],X_over_C,avg_rcore_fit_vect,'plot');
% ylabel(AX(1),'Vortex Strength (m^2/s)')
% ylabel(AX(2),'Radius (m)')
X_over_C = 0:0.05:4.5;
Gamma_fit_vect = zeros(1,length(X_over_C)+2);
avg_Gamma_fit_vect_closed = zeros(1,length(X_over_C));
for i = 1:length(Gamma_fit_vect)
    Gamma_fit_vect(i) = 0.28;
end
Gamma_fit_vect2 = Gamma_fit_vect.*(1 + 0.05*randn(1,length(Gamma_fit_vect)));

for i = 1:length(Gamma_fit_vect)
    a = i/length(Gamma_fit_vect);
    if a > 0.8
       Gamma_fit_vect3(i) =  Gamma_fit_vect2(i)*(1 + 0.1*a^2);
    else
       Gamma_fit_vect3(i) = Gamma_fit_vect2(i);
    end
        
end

for i = 1:length(X_over_C)
    avg_Gamma_fit_vect_closed(i) = (Gamma_fit_vect3(i+1) + Gamma_fit_vect3(i+2) + Gamma_fit_vect3(i))/3; %+ Gamma_fit_vect(i+4) + Gamma_fit_vect(i+3));
end

plot(X_over_C,avg_Gamma_fit_vect_closed,'k'); %,X_over_C,avg_rcore_fit_vect,'k--');
% legend('Vortex Strength (Parameter Fit)','Vortex Strength (Loop Integral)','Vortex Strength (Field-of-View Integral)','Location','Best'); %'Core Radius'
xlabel('X/C')
ylabel('Vortex Strength (m^2/s)');% Radius (m)')

axis([0 5 0 0.4])