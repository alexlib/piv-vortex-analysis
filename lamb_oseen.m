function F = lamb_oseen(coeffs,r)

twopi = 2*pi;
r(r==0) = inf;
F = (coeffs(1)./(twopi.*(r))).*(1-exp((-r.*r)./(coeffs(2)*coeffs(2))));

end
