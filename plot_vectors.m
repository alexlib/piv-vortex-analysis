function [] = plot_vectors(mean_gamma_grid,gamma_set,length_unique_x,length_unique_y,window_size,loc_x,loc_y,x,y,u,v,Px,Py,unique_x,unique_y,P_pos_temp,clearance,gamma_max,name)

%REARRANGE GAMMA DATA SO THAT IT CAN BE PLOTTED
gamma_line = gamma_set(P_pos_temp,:);
gamma_dist = zeros(length_unique_x,length_unique_y);
counter1 = 0;

for ix = 1:window_size
    for iy = 1:window_size
        counter1 = iy+window_size*(ix-1);     
        gamma_dist(iy + loc_y - clearance - 1,ix + loc_x - clearance) = gamma_line(counter1);
    end
end

%PRODUCE PLOT OF DATA
quiver(x,y,u,v);
xlabel('X(m)');
ylabel('Y(m)');
hold on;
plot(Px,Py,'kx'); %Plot actual vortex location
hold on;
plot(Px,Py,'ko'); %Plot detected vortex location
% V = 0.6*gamma_max:0.1*gamma_max:gamma_max; %Set gamma values for contour plotting

% contour(unique_x,unique_y,mean_gamma_grid); %,V);

title(name)


hold on;
line('XData', [-0.004 0.1603], 'YData', [.003 0.003], 'LineStyle', '-', 'LineWidth', 2, 'Color','k');
hold on;
line('XData', [-0.004 -0.004], 'YData', [.003 0.1647], 'LineStyle', '-', 'LineWidth', 2, 'Color','k');
hold on;
line('XData', [0.002 0.1603], 'YData', [0.1 0.1], 'LineStyle', '-', 'LineWidth', 2, 'Color','k');
hold on;
line('XData', [0.002 0.002], 'YData', [.003 .1647], 'LineStyle', '-', 'LineWidth', 2, 'Color','k');
hold on;
line('XData', [0.002 0.1603], 'YData', [0.065 0.065], 'LineStyle', '-', 'LineWidth', 1, 'Color','k');
hold on;
line('XData', [0.002 0.1603], 'YData', [0.04 0.04], 'LineStyle', '-', 'LineWidth', 1, 'Color','k');


hold off

end

