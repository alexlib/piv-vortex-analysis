function[Gamma_FOV,ndx,dx,ndy,dy] = integral_fov(xmax,xmin,length_unique_x,ymax,ymin,length_unique_y,xyuv,N)

%DEFINE UNIVERSAL CONSTANTS
ndx = length_unique_x - 1;
ndy = length_unique_y - 1;
dx = (xmax - xmin)/ndx;
dy = (ymax - ymin)/ndy;

%INTEGRAL FOR (X, YMIN)
u_ymin = zeros(ndx,1);
for i = 1:ndx
    u_ymin(i,1) = 0.5*(xyuv((i)*length_unique_y+1,3) + xyuv((i-1)*length_unique_y+1,3));
end
int_ymin = sum(u_ymin)*dx;

%INTEGRAL FOR (X, YMAX)
u_ymax = zeros(ndx,1);
for i = 1:ndx
    u_ymax(i,1) = -0.5*(xyuv((i+1)*length_unique_y-1,3) + xyuv(i*length_unique_y-1,3));
end
int_ymax = sum(u_ymax)*dx;

%INTEGRAL FOR (XMIN, Y)
v_xmin = zeros(ndy,1);
for i = 1:ndy
    v_xmin(i,1) = -0.5*(xyuv(i+1,4) + xyuv(i,4));
end
int_xmin = sum(v_xmin)*dy;

%INTEGRAL FOR (XMAX, Y)
v_xmax = zeros(ndy,1);
for i = 1:ndy
    v_xmax(i,1) = 0.5*(xyuv(N-length_unique_y+i+1,4) + xyuv(N-length_unique_y+i,4));
end
int_xmax = sum(v_xmax)*dy;

%SUM COMPONENTS
Gamma_FOV = int_xmax + int_xmin + int_ymax + int_ymin;

end