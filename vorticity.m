function [omega_max,vort_centre_row,vort_centre_column,vort_centre_x,vort_centre_y] = vorticity(ndy,ndx,umesh,vmesh,dx,dy,xmin,ymin,ymax)

%CALCULATE VORTICITY FOR EACH POINT ON GRID
dv_x = zeros(ndy,ndx);
du_y = zeros(ndy,ndx);
omega = zeros(ndy,ndx);
for ix = 1:ndx
    for iy = 1:ndy
        du_y(iy,ix) = umesh(iy+1,ix) - umesh(iy,ix);
        dv_x(iy,ix) = vmesh(iy,ix+1) - vmesh(iy,ix);
        omega(iy,ix) = dv_x(iy,ix)/dx - du_y(iy,ix)/dy;
    end
end

%DETERMINE LOCATION AND VALUE OF MAXIMUM VORTICITY
[omega_max, vort_centre_temp] = max(abs(omega(:)));
[vort_centre_row,vort_centre_column] = ind2sub(size(omega),vort_centre_temp);

vort_centre_x = xmin + dx*(vort_centre_column - 1);
vort_centre_y = ymax - dy*(vort_centre_row - 1);
end

