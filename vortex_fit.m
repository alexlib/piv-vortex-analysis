function [v_theta,Gamma_fit,rcore_fit] = vortex_fit(Gamma0,rcore0,xyuv,Px,Py,length_xmesh,length_ymesh,umesh,vmesh,xmesh,ymesh)

%REARRANGE DATA INTO GRID
% [xmesh,ymesh] = meshgrid(xmin:(xmax-xmin)/(length_unique_x-1):xmax,ymax:(ymin-ymax)/(length_unique_y-1):ymin);
% length_xmesh = length(xmesh);
% length_ymesh = length(ymesh);
% umesh = zeros(length_ymesh,length_xmesh);
% vmesh = zeros(length_ymesh,length_xmesh);
rsq = zeros(length_ymesh,length_xmesh);
r = zeros(length_ymesh,length_xmesh);
theta = zeros(length_ymesh,length_xmesh);

for ix = 1:length_xmesh
    for iy = 1:length_ymesh
        counter = iy+length_ymesh*(ix-1);     
%         umesh(length_ymesh+1-iy,ix) = xyuv(counter,3);
%         vmesh(length_ymesh+1-iy,ix) = xyuv(counter,4);
        rsq(iy,ix) = (xmesh(iy,ix)-Px)*(xmesh(iy,ix)-Px)+(ymesh(iy,ix)-Py)*(ymesh(iy,ix)-Py);
        r(iy,ix) = sqrt(rsq(iy,ix));
        theta(iy,ix) = atan2(ymesh(iy,ix)-Py,xmesh(iy,ix)-Px);
    end
end

%DETERMINE SIGN OF V_THETA
% v_theta_sign = zeros(length_xmesh,length_ymesh);
% for ix = 1:length_xmesh
%     for iy = 1:length_ymesh
%         if theta(iy,ix) > -0.5*pi && theta(iy,ix) < 0.5*pi
%         v_theta_sign(iy,ix) = sign(vmesh(iy,ix));
%         elseif theta(iy,ix) == 0.5*pi
%         v_theta_sign(iy,ix) = -sign(umesh(iy,ix));    
%         elseif theta(iy,ix) == -0.5*pi
%         v_theta_sign(iy,ix) = sign(umesh(iy,ix));      
%         elseif theta(iy,ix) == pi
%         v_theta_sign(iy,ix) = -sign(vmesh(iy,ix));
%         else
%         v_theta_sign(iy,ix) = -sign(vmesh(iy,ix));
%         end
%     end
% end

%DETERMINE TANGENTIAL VELOCITIES
v_thetasq = umesh.*umesh + vmesh.*vmesh;
v_theta = sqrt(v_thetasq);
% v_theta = sqrt(v_thetasq).*v_theta_sign;

%CONVERT FROM MATRICES INTO VECTORS
r_vect = zeros(1,length_xmesh*length_xmesh);
v_theta_vect = zeros(1,length_xmesh*length_xmesh);
for ix = 1:length_xmesh
    for iy = 1:length_ymesh
        counter = iy+length_xmesh*(ix-1);
        r_vect(1,counter) = r(iy,ix);
        v_theta_vect(1,counter) = v_theta(iy,ix);
    end
end

%PERFORM LEAST SQUARES CALCULATIONS TO DETERMINE GAMMA AND CORE RADIUS
coeffs0 = [Gamma0; rcore0]; %Initial guess
% OPTIONS = optimset('Algorithm','levenberg-marquardt',@lsqcurvefit);
coeffs = lsqcurvefit(@lamb_oseen,coeffs0,r_vect,v_theta_vect);
Gamma_fit = coeffs(1);
rcore_fit = coeffs(2);


