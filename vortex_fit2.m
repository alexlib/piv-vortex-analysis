function [v_theta,Gamma_fit,rcore_fit] = vortex_fit2(Gamma0,rcore0,Px,Py,umesh,vmesh,xmesh,ymesh,xpos,ypos)

box_size = 15;
plus_min = (box_size - 1)/2;

ubox = zeros(box_size);
vbox = zeros(box_size);
rsq = zeros(box_size);
r = zeros(box_size);
% theta = zeros(box_size);

for ix = 1:box_size
    for iy = 1:box_size
%         counter = iy+box_size*(ix-1);
        cy = ypos - plus_min + iy - 1;
        cx = ix + xpos - plus_min - 1;
        ubox(iy,ix) = umesh(cy,cx);
        vbox(iy,ix) = vmesh(cy,cx);
       
        rsq(iy,ix) = (xmesh(cy,cx)-Px)*(xmesh(cy,cx)-Px)+(ymesh(cy,cx)-Py)*(ymesh(cy,cx)-Py);
        r(iy,ix) = sqrt(rsq(iy,ix));
%         theta(iy,ix) = atan2(ymesh(cy,cx)-Py,xmesh(cy,cx)-Px);
    end
end

%DETERMINE TANGENTIAL VELOCITIES
v_thetasq = ubox.*ubox + vbox.*vbox;
v_theta = sqrt(v_thetasq);

%CONVERT FROM MATRICES INTO VECTORS
box_sizesq = box_size*box_size;
r_vect = zeros(1,box_sizesq);
v_theta_vect = zeros(1,box_sizesq);
for ix = 1:box_size
    for iy = 1:box_size
        counter = iy+box_size*(ix-1);
        r_vect(1,counter) = r(iy,ix);
        v_theta_vect(1,counter) = v_theta(iy,ix);
    end
end

%PERFORM LEAST SQUARES CALCULATIONS TO DETERMINE GAMMA AND CORE RADIUS
coeffs0 = [Gamma0; rcore0]; %Initial guess
lb = [0; 0];
ub = [1; 0.5];
OPTIONS = optimset('Algorithm','trust-region-reflective');
coeffs = lsqcurvefit(@lamb_oseen,coeffs0,r_vect,v_theta_vect,lb,ub,OPTIONS);
Gamma_fit = coeffs(1);
rcore_fit = coeffs(2);

end

