function [] = vortex_generator(vortex_strength,grid_density,core_radius,x_coord,y_coord)

%***************************************************
%                  GENERATE VORTEX
%***************************************************

%GENERATE GRID
xmin = -0.2;
xmax = 0.2;
ymin = -0.2;
ymax = 0.2;
[x,y] = meshgrid(xmin:((xmax-xmin)/grid_density):xmax,ymin:((ymax-ymin)/grid_density):ymax);
length_x = length(x);
length_y = length(y);

%GENERATE VELOCITIES, TAKING CENTRE AS (0,0)
r = zeros(length_x,length_y);
rsq = zeros(length_x,length_y);
theta = zeros(length_x,length_y);
v_theta = zeros(length_x,length_y);
u = zeros(length_x,length_y);
v = zeros(length_x,length_y);

%BRING TRIVIAL CALCULATIONS OUTSIDE OF LOOP FOR EFFICIENCY
twopi = 2*pi;
% four_vt = (core_radius/1.121)*(core_radius/1.121);

%GENERATE THE VORTEX VELOCITY VECTORS
for ix = 1:length_x
    for iy = 1:length_y
        rsq(ix,iy) = (x(ix,iy)-x_coord)*(x(ix,iy)-x_coord)+(y(ix,iy)-y_coord)*(y(ix,iy)-y_coord);
        r(ix,iy) = sqrt(rsq(ix,iy));
        theta(ix,iy) = atan2(y(ix,iy)-y_coord,x(ix,iy)-x_coord);
        
        %IDEAL VORTEX
        %v_theta(ix,iy) = gamma/(twopi*r(ix,iy));
        
        %LAMB-OSEEN VORTEX
        v_theta(ix,iy) = (vortex_strength/(twopi*r(ix,iy)))*(1-exp(-rsq(ix,iy)/(core_radius*core_radius)));
        
        %CALCULATE U AND V VELOCITY COMPONENTS
        if -0.5*pi <= theta(ix,iy) <= 0.5*pi
        u(ix,iy) = -v_theta(ix,iy)*sin(theta(ix,iy));
        v(ix,iy) = v_theta(ix,iy)*cos(theta(ix,iy));
        else
        u(ix,iy) = v_theta(ix,iy)*sin(theta(ix,iy));
        v(ix,iy) = v_theta(ix,iy)*cos(theta(ix,iy));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ADD GAUSSIAN NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%KEEP TRACK OF ORIGINAL VELOCITIES FOR TESTING PURPOSES
% u_clean = u;
% v_clean = v;

%ADD GAUSSIAN NOISE
u = u.*(1+0.2*randn(length_x,length_y));
v = v.*(1+0.2*randn(length_x,length_y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     END OF ADD GAUSSIAN NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT RESULTS TO A TEXT FILE
output_array = zeros(4,length_x*length_y);

for ix = 1:length_x
    for iy = 1:length_y
        counter = iy+length_x*(ix-1);
        output_array(1,counter) = x(ix,iy);
        output_array(2,counter) = y(ix,iy);
        output_array(3,counter) = u(ix,iy);
        output_array(4,counter) = v(ix,iy);
    end
end

vortex = fopen('vortex.txt','w');
fprintf(vortex, 'X Y U V a b c d e f g h \r\n');
fprintf(vortex, '%f %f %f %f\r\n',output_array);
fclose(vortex);

%DRAW VELOCITY VECTORS ON THE SCREEN
quiver(x,y,u,v)
xlabel('X(m)')
ylabel('Y(m)')

%***************************************************
%               END OF GENERATE VORTEX
%***************************************************

end

