%*************************************************************************
%                    START OF MAIN PROGRAM
%*************************************************************************

%SIGN CONVENTION - X GOES LEFT TO RIGHT (-5 TO 5), Y GOES DOWN TO UP (-5 TO 5)
%UNITS - VELOCITIES ARE IN METRES PER SECOND. DISTANCES ARE IN M.

%CLEAR ALL VARIABLES
clear

%*************************************************************************
%           ACTUAL DATA: READ DATA IN FROM DATA FILE 'vortex.txt'
%*************************************************************************

% DEFINE NAME OF DATA SET
name = '2msSealed-F0636-500DP-300dt';
no_of_files = 3;
grid_density_real = 64; %64
grid_density_real_sq = grid_density_real*grid_density_real;

% DEFINE RELEVANT FRAMES TO ANALYSE:
if     strcmp(name,'2msSealed-50O370-500DP-300dt')
    startframes = [39 37 34 31 45];
elseif strcmp(name,'2msSealed-FO545-500DP-300dt')
    startframes = [34 41 31 45 36];
elseif strcmp(name,'2msSealed-FC195-500DP-300dt')
    startframes = [32 29 32 49 41];
elseif strcmp(name,'2msSealed-250283-500DP-300dt')
    startframes = [46 60 51 50 37];
elseif strcmp(name,'2msSealed-750458-500DP-300dt')
    startframes = [53 47 53 40 79];
elseif strcmp(name,'2msSealed-F0636-500DP-300dt')
    startframes = [45 56 44];
elseif strcmp(name,'2msSealed-FO720-500DP-300dt')
    startframes = [50 44 83];
elseif strcmp(name,'3msSealed-FC-500DP-300dt')
    startframes = [40]; %40
elseif strcmp(name,'FC195-1000SP')
    startframes = [67 73 64 78 68];
elseif strcmp(name,'FCNoCam-1000SP')
    startframes = [63 112 61 58 67];
elseif strcmp(name,'FC195-1000SP')
    startframes = [67 73 64 78 68];
else
    disp('Error')
end

length_startframes = length(startframes);

savefile = strcat(name, '.txt');
results = fopen(savefile,'w');

%LOOP SEVERAL TIMES TO RUN DOWNSTREAM
zmax = 2; %98
startframes0 = startframes + 0;

for z = 1:zmax

startframes = startframes0 + z-1 -0   ;  %Add frames downstream

% DEFINE VARIABLES
xyuv_noise = zeros(grid_density_real_sq,no_of_files*4);
xyuv_noise_unsorted = zeros(grid_density_real_sq,no_of_files*4);
filestart = cell(1,no_of_files);

% READ IN EACH DATA SET
for k = 1:no_of_files
%     READ IN EACH DATA SET
    if startframes(k) >= 100
        filestart{k} = 'B00';
    else
        filestart{k} = 'B000';
    end
    cd(strcat(name,'-',num2str(k),'_PIV_Vec_PreProc_MP(32x32_50ov)_PostProc=unknown'));
    input_filename = strcat(filestart{k},num2str(startframes(k)),'.txt');
    vortex = fopen(input_filename);
    titles = fscanf(vortex, '%s %s %s %s %s %s %s %s %s %s %s %s', [12 1]);  
    k_low = (k-1)*4 + 1  ;
    k_high = k*4 ;
    xyuv_noise_unsorted(:,k_low:k_high) = fscanf(vortex, '%f %f %f %f', [4 inf])';
    xyuv_noise(:,k_low:k_high) = sortrows(xyuv_noise_unsorted(:,k_low:k_high));
    fclose(vortex);
    cd ..
end



%COMBINE FILES TO REDUCE NOISE
no_zeros = zeros(grid_density_real_sq,1);
xyuv = zeros(grid_density_real_sq,4);
total_different_velocity_components = no_of_files*2;

%%%%%%%%%%%%%%%%%% NOTE: THIS LOOP IS AFFECTED BY K %%%%%%%%%%%%%%
for i = 1:grid_density_real_sq
    no_zeros(i,1) = sum(xyuv_noise(i,:) == 0);
    xyuv(i,3:4) = (  xyuv_noise(i,3:4) + xyuv_noise(i,7:8) + xyuv_noise(i,11:12)  )/((total_different_velocity_components - no_zeros(i,1))*0.5);
     %xyuv_noise(i,3:4) + xyuv_noise(i,7:8) + xyuv_noise(i,11:12)  + xyuv_noise(i,15:16) + xyuv_noise(i,19:20)
     
     if abs(xyuv(i,3)) > 2.2
         xyuv(i,3) = 0;
         xyuv(i,4) = 0;
     end
     if abs(xyuv(i,4)) > 2.2
         xyuv(i,4) = 0;
         xyuv(i,3) = 0;
     end
end
xyuv(:,1:2) =  xyuv_noise(:,1:2)/1000; %Converting to metres

%*************************************************************************
%           SYNTHETIC DATA: READ DATA IN FROM DATA FILE 'vortex.txt'
%*************************************************************************
% 
% % GENERATE LAMB-OSEEN VORTEX DATASET:
% % vortex_generator(vortex strength,grid_density,core_radius,x_coord,y_coord)
% % Works up to 5 steps from xmax,ymax
% vortex_generator(-0.4,60,0.04,0.02,-0.02);
% 
% %OPEN THE FILE AND READ IN THE DATA
% vortex = fopen('vortex.txt');
% titles = fscanf(vortex, '%s %s %s %s %s %s %s %s %s %s %s %s', [12 1]); % 
% xyuv_temp = fscanf(vortex, '%f %f %f %f', [4 inf])';
% xyuv = sortrows(xyuv_temp);
% fclose(vortex);
% name = 'results';

%*************************************************************************
%                           PROCESS DATA
%*************************************************************************

%CONVERT NaN VALUES TO ZERO (ASSUMING THAT THE ONLY NAN VALUES ARE THE
%VELOCITIES AT THE CENTRE OF THE GRID)
xyuv(isnan(xyuv)) = 0;
xyuv(xyuv == 0) = randn/100000; %Necessary to avoid the areas with no velocity giving gamma_max = 1.0

%CONVERT TEXT FILE DATA INTO VECTORS
x = xyuv(:,1)'; y = xyuv(:,2)'; u = xyuv(:,3)'; v = xyuv(:,4)'; M = xyuv(:,1:2)'; P = xyuv(:,1:2)'; U_M = xyuv(:,3:4)';

%REQUIRED CONSTANTS
unique_x = unique(x); unique_y = unique(y);
length_unique_x = length(unique_x); length_unique_y = length(unique_y);
N = length(xyuv); Nsq = N*N;
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);

%*************************************************************************
%                          REARRANGE DATA
%*************************************************************************

%REARRANGE DATA INTO GRID
[xmesh,ymesh] = meshgrid(xmin:(xmax-xmin)/(length_unique_x-1):xmax,ymax:(ymin-ymax)/(length_unique_y-1):ymin);
length_xmesh = length(xmesh); length_ymesh = length(ymesh);
umesh = zeros(length_ymesh,length_xmesh);
vmesh = zeros(length_ymesh,length_xmesh);

for ix = 1:length_xmesh
    for iy = 1:length_ymesh
        counter = iy+length_ymesh*(ix-1);     
        umesh(length_ymesh+1-iy,ix) = xyuv(counter,3);
        vmesh(length_ymesh+1-iy,ix) = xyuv(counter,4);
    end
end

%*************************************************************************
%                CALCULATE GAMMA FOR EACH POINT ON GRID
%*************************************************************************

% CALL GAMMA FUNCTION TO FIND VALUE OF GAMMA FOR EACH LOCATION, P
window_size = 11;
[gamma,Np,window_size_sq,clearance] = gamma_function(window_size,length_unique_x,length_unique_y,M,P,U_M);

%*************************************************************************
%   DETERMINE LOCATION OF VORTEX BY FINDING WHICH LOCATION OF P GIVES THE
%   LARGEST AVERAGE GAMMA VALUE
%*************************************************************************

%STORE ALL VALUES OF GAMMA FOR EACH P IN A SEPARATE ROW
gamma_set = zeros(Np, window_size_sq);
mean_gamma = zeros(Np,1);
for i = 1:Np
    i_low = (i-1)*window_size_sq + 1;
    i_high = i*window_size_sq;
    gamma_set(i,:) = gamma(i_low:i_high);
    mean_gamma(i) = mean(gamma_set(i,:));
end

%DETERMINE LOCATION OF VORTEX (Px,Py)
[gamma_max, P_pos_temp] = max(abs(mean_gamma));
loc_x = clearance + ceil(P_pos_temp/(length_unique_y - 2*clearance)); %Rounding up
P_pos = P_pos_temp + clearance*length_unique_y + 2*clearance*(loc_x - clearance) - clearance; %Locate position of P in XYUV, factoring in for areas where P was not tested for
Px = x(P_pos);
loc_y = mod(P_pos,length_unique_x);
Py = y(loc_y);
xpos = loc_x;
ypos = length_unique_y - loc_y + 1;

mean_gamma_grid = zeros(length_unique_y,length_unique_x);
for iPx = 1:(length_unique_x - 2*clearance)
    for iPy = 1:(length_unique_y - 2*clearance)
        counter = iPy + (length_unique_y - 2*clearance)*(iPx-1);
        mean_gamma_grid(clearance+iPy,clearance+iPx) = mean_gamma(counter);
    end
end

%*************************************************************************
%   DETERMINE CIRCULATION OF VORTEX BY INTEGRATING AROUND FIELD OF VIEW
%*************************************************************************

[Gamma_FOV,ndx,dx,ndy,dy] = integral_fov(xmax,xmin,length_unique_x,ymax,ymin,length_unique_y,xyuv,N);

%*************************************************************************
%   DETERMINE CIRCULATION OF VORTEX USING VORTEX PARAMETER FITTING
%*************************************************************************

if abs(Px) <= 0.1 && abs(Py) <= 0.1
    [v_theta,Gamma_fit,rcore_fit] = vortex_fit3(-0.3,0.03,Px,Py,umesh,vmesh,xmesh,ymesh,xpos,ypos);
else
    Gamma_fit = Gamma_fit;
end

% if rcore_fit >= 0.05
%     break
% end

if rcore_fit >= 0.05
    rcore_fit = 0.05;
end

%*************************************************************************
%   VERIFY LOCATION OF CENTRE OF VORTEX BY DETERMINING VORTICITY MAXIMUM
%*************************************************************************

[omega_max,vort_centre_row,vort_centre_column,vort_centre_x,vort_centre_y] = vorticity(ndy,ndx,umesh,vmesh,dx,dy,xmin,ymin,ymax);

%*************************************************************************
%   DETERMINE CIRCULATION OF VORTEX BY INTEGRATING AROUND BOX OF SIDE
%   LENGTH 2.5*RCORE_FIT
%*************************************************************************

if abs(Px) <= 0.1 && abs(Py) <= 0.1
    [Gamma_box,ndx_box,ndy_box] = integral_box(rcore_fit,dx,dy,xpos,ypos,umesh,vmesh,2);
else
    Gamma_box = Gamma_fit;
end
    

% OUTPUT DATA TO SCREEN

% sprintf('VORTEX CENTRE: Gamma Method = (%d,%d)m; Vorticity Method = (%d,%d)m',Px,Py,vort_centre_x,vort_centre_y)
% sprintf('VORTEX CORE RADIUS: Parameter Fit = %dm',(rcore_fit))
% sprintf('VORTEX STRENGTH: Parameter Fit = %dm^2/s; Integral FOV = %dm^2/s, Integral Box = %dm^2/s',Gamma_fit,Gamma_FOV,Gamma_box)

%WRITE DATA TO FILE

% savefile = strcat(name, '.txt');
% results = fopen(savefile,'w');
fprintf(results,'FRAME NUMBER = %d \r\n',z);
fprintf(results,'VORTEX CENTRE: Gamma Method = (%d,%d); Vorticity Method = (%d,%d)\r\n',Px,Py,vort_centre_x,vort_centre_y);
fprintf(results,'VORTEX CORE RADIUS: Parameter Fit = %dm \r\n',abs(rcore_fit));
fprintf(results,'VORTEX STRENGTH: Parameter Fit = %dm^2/s; Integral FOV (Ineffective due to noise) = %dm^2/s, Integral Box = %dm^2/s \r\n \r\n',Gamma_fit,Gamma_FOV,Gamma_box);
% fclose(results);

%*************************************************************************
%                        DATA ANALYSIS
%*************************************************************************

figure;
plot_vectors(mean_gamma_grid,gamma_set,length_unique_x,length_unique_y,window_size,loc_x,loc_y,x,y,u,v,Px,Py,unique_x,unique_y,P_pos_temp,clearance,gamma_max,name);

%PLOT VARIATION OF PARAMETERS WITH DISTANCE DOWNSTREAM
Gamma_fit_vect(z) = abs(Gamma_fit);
Gamma_box_vect(z) = abs(Gamma_box);
rcore_fit_vect(z) = rcore_fit;
Px_vect(z) = Px;
Py_vect(z) = Py;
t(z) = (z-1)*0.002;

end

if rcore_fit >= 0.05 && z > 1;
    Gamma_fit_vect(z) = Gamma_fit_vect(z-1);
    Gamma_box_vect(z) = Gamma_box_vect(z-1);
    rcore_fit_vect(z) = rcore_fit_vect(z-1);
    Px_vect(z) = Px;
    Py_vect(z) = Py;
    t(z) = t(z-1)+0.002;
end

X = t.*2;
X_over_C = X/0.08; %Non-dimensionalised using chord length of flap

% % PLOT CIRCULATION AND CORE RADIUS AGAINST TIME
% figure;
% plot(X_over_C,Gamma_fit_vect,'k',X_over_C,Gamma_box_vect,'b',X_over_C,rcore_fit_vect,'k--');
% legend('Vortex Strength (Parameter Fit)','Vortex Strength (Integral)','Core Radius','Location','Best')
% xlabel('X/C')
% ylabel('Vortex Strength (m^2/s); Radius (m)')
% % ylim(
% title(name)

%PLOT VARIATION IN CENTRE OF VORTEX
% gradient_x = (Px_vect(z) - Px_vect(1))/X(z);
% gradient_y = (Py_vect(z) - Py_vect(1))/X(z);
% theta_x = 180*(atan2(Px_vect(z) - Px_vect(1),X(z)))/pi;
% theta_y = 180*(atan2(Py_vect(z) - Py_vect(1),X(z)))/pi;
% str1(1) = {['Theta_x = ',num2str(theta_x)]};
% str1(2) = {['Theta_y = ',num2str(theta_y)]};
% plot(X,Px_vect,'k',X,Py_vect,'k--');
% legend('Px','Py')
% xlabel('Distance (m)')
% ylabel('Location (m)')
% text(50,20,str1)
% title(name)

%*************************************************
%             END OF MAIN PROGRAM
%*************************************************

