function [v_theta,Gamma_fit,rcore_fit] = vortex_fit4(Gamma0,rcore0,Px,Py,umesh,vmesh,xmesh,ymesh,xpos,ypos,nrcore_fit,z)




box_size = ceil(1.25*nrcore_fit);

plus_min = ceil(box_size/2);

ubox = zeros(box_size);
vbox = zeros(box_size);
rsq = zeros(box_size);
r = zeros(box_size);
theta = zeros(box_size);

for ix = 1:box_size
    for iy = 1:box_size
        cy = ypos - plus_min + iy - 1;
        cx = ix + xpos - plus_min - 1;
        ubox(iy,ix) = umesh(cy,cx);
        vbox(iy,ix) = vmesh(cy,cx);
        xbox(iy,ix) = xmesh(cy,cx);
        ybox(iy,ix) = ymesh(cy,cx);
        rsq(iy,ix) = (xmesh(cy,cx)-Px)*(xmesh(cy,cx)-Px)+(ymesh(cy,cx)-Py)*(ymesh(cy,cx)-Py);
        theta(iy,ix) = atan2(ymesh(cy,cx)-Py,xmesh(cy,cx)-Px); %OK
    end
end

r = sqrt(rsq);

%DETERMINE SIGN OF V_THETA
v_theta_sign = zeros(box_size);

for ix = 1:box_size
    for iy = 1:box_size
        if theta(iy,ix) > -0.5*pi && theta(iy,ix) < 0.5*pi
        v_theta_sign(iy,ix) = sign(vbox(iy,ix));
        elseif theta(iy,ix) == 0.5*pi
        v_theta_sign(iy,ix) = -sign(ubox(iy,ix));    
        elseif theta(iy,ix) == -0.5*pi
        v_theta_sign(iy,ix) = sign(ubox(iy,ix));      
        elseif theta(iy,ix) == pi
        v_theta_sign(iy,ix) = -sign(vbox(iy,ix));
        else
        v_theta_sign(iy,ix) = -sign(vbox(iy,ix));
        end
    end
end

%DETERMINE TANGENTIAL VELOCITIES
v_thetasq = ubox.*ubox + vbox.*vbox;
v_theta = sqrt(v_thetasq).*v_theta_sign;

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
lb = [-0.5; 0];
ub = [0; 0.05];
OPTIONS = optimset('Algorithm','trust-region-reflective');
coeffs = lsqcurvefit(@lamb_oseen,coeffs0,r_vect,v_theta_vect,lb,ub,OPTIONS);
Gamma_fit = coeffs(1);
rcore_fit = coeffs(2);

end

