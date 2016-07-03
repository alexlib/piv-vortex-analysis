function [Gamma_box,ndx_box,ndy_box] = integral_box(rcore_fit,dx,dy,xpos,ypos,umesh,vmesh,size)


%NUMBER OF DIVISIONS TO LOOK OUT
% ndx_box = floor(abs(rcore_fit)/dx);
% ndy_box = floor(abs(rcore_fit)/dy);
% shift_y = floor((ndy_box)/2);
% shift_x = floor((ndx_box1)/2);
shift_x = ceil(abs(size*0.5*1.121*rcore_fit)/dx);
shift_y = ceil(abs(size*0.5*1.121*rcore_fit)/dy);
% shift_x = 5;
% shift_y = 5;
ndx_box = 2*shift_x + 1;
ndy_box = 2*shift_y + 1;

% INTEGRAL FOR (X, YLOW)
u_ylow = zeros(ndx_box,1);
for i = 1:ndx_box
    u_ylow(i,1) = 0.5*(umesh(ypos+shift_y,xpos-shift_x+1+i) + umesh(ypos+shift_y,xpos-shift_x+i));
end
int_ylow = sum(u_ylow)*dx;

% INTEGRAL FOR (X, YHIGH)
u_yhigh = zeros(ndx_box,1);
for i = 1:ndx_box
    u_yhigh(i,1) = -0.5*(umesh(ypos-shift_y,xpos-shift_x+1+i) + umesh(ypos-shift_y,xpos-shift_x+i));
end
int_yhigh = sum(u_yhigh)*dx;

% INTEGRAL FOR (XLOW, Y)
v_xlow = zeros(ndy_box,1);
for i = 1:ndy_box
    v_xlow(i,1) = -0.5*(vmesh(ypos-shift_y+1+i,xpos-shift_x) + vmesh(ypos-shift_y+i,xpos-shift_x));
end
int_xlow = sum(v_xlow)*dy;

% INTEGRAL FOR (XHIGH, Y)
v_xhigh = zeros(ndy_box,1);
for i = 1:ndy_box
    v_xhigh(i,1) = 0.5*(vmesh(ypos-shift_y+1+i,xpos+shift_x) + vmesh(ypos-shift_y+i,xpos+shift_x));
end
int_xhigh = sum(v_xhigh)*dy;

Gamma_box = int_ylow + int_yhigh + int_xlow + int_xhigh;

end

