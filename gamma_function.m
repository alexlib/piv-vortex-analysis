function [gamma,Np,window_size_sq,clearance] = gamma_function(window_size,length_unique_x,length_unique_y,M,P,U_M)

%SET CONVECTION VELOCITY
U_P = [0; 0]; %Be careful when choosing U_P so as not to make it too large

%DEFINE WINDOW SIZE
% window_size = 11;
window_size_sq = window_size*window_size;
clearance = (window_size - 1)/2; %Number of points in x,y directions around considering point, P
clearancepl1 = clearance + 1;
Np = (length_unique_x - 2*clearance)*(length_unique_y - 2*clearance); %Total number of locations tested for P

%CALCULATE GAMMA, ITERATING M WITHIN P
% PM = zeros(3,Nsq);
% U_M_minus_U_P = zeros(3,Nsq);
% cross_product_dot_z = zeros(1,Nsq);
% dot_product = zeros(1,Nsq);
% gamma = zeros(1,Nsq);

PM = zeros(3,Np*window_size_sq);
U_M_minus_U_P = zeros(3,Np*window_size_sq);
cross_product_dot_z = zeros(1,Np*window_size_sq);
dot_product = zeros(1,Np*window_size_sq);
gamma = zeros(1,Np*window_size_sq);
counter = 0;

for iPx = clearancepl1:(length_unique_x - clearance)
    for iPy = clearancepl1:(length_unique_y - clearance)
        
        %DETERMINE CORRECT POSITION IN XYUV MATRIX
        iP = (iPx-1)*length_unique_y + iPy;
        
        for iMx = (iPx - clearance):(iPx + clearance)
            for iMy = (iPy - clearance):(iPy + clearance)
                iM = (iMx-1)*length_unique_y + iMy;
%                 counter = iM+N*(iP-1);
                counter = counter + 1;
                PM(1:2,counter) = M(:,iM) - P(:,iP);
                U_M_minus_U_P(1:2,counter) = U_M(:,iM); %Zero convection velocity so no need for - U_P;

                %Cross product will only fill 3rd row due to 2D flow so write as multiplication to speed up calculations
                cross_product_dot_z(counter) = PM(1,counter)*U_M_minus_U_P(2,counter) - PM(2,counter)*U_M_minus_U_P(1,counter);

                %MATLAB DOT function involves unnecessary checks - rewriting the base code speeds up calculations
        %         dot_product(counter) = dot(norm(PM(:,counter)),norm(U_M_minus_U_P(:,counter)));
                dot_product(counter) = sum(conj(norm(PM(:,counter))).*norm(U_M_minus_U_P(:,counter)));

                gamma(counter) = (cross_product_dot_z(counter)/dot_product(counter));
            end
        end  
    end
end

%CONVERT NaN VALUES (NaN values arise from where PM = 0)
gamma(isnan(gamma)) = 1;

end

