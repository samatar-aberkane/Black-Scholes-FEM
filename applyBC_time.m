function [L_bc, b_bc] = applyBC_time(L, b, dFreedom, Smax, K, r, tau)
    m = max(max(dFreedom));
    L_bc = L;
    b_bc = b;

    u = zeros(m,1);
    u(1)   = 0;
    u(m)   = Smax - K * exp(-r * tau);  % time-dependent right BC

    % left DOF index (usually 1)
    cL = 1;
    b_bc = b_bc - L_bc*u;
    b_bc(cL) = u(cL);
    L_bc(cL,:) = 0;
    L_bc(:,cL) = 0;
    L_bc(cL,cL) = 1;

    % right DOF index (usually m)
    cR = m;
    b_bc(cR) = u(cR);
    L_bc(cR,:) = 0;
    L_bc(:,cR) = 0;
    L_bc(cR,cR) = 1;
end
