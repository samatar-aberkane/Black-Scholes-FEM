function [kQ,cQ,bQ,fQ]=evals(xq)
E = 5;
rho_g = 10;   % = rho * g
x_0 = 1;
l = 5;

A = (H(xq) - H(xq - x_0)) + (xq*(2*l - xq)/l.^2) * H(xq - x_0);

kQ = E*A;       % coefficient of (u')' term in strong form
cQ = 0;         % no convection-like term
bQ = 0;         % no reaction term
fQ = -rho_g;    % body force
end

function [out]=H(x)
if x < 0
    out=0;
else
    out=1;
end
end