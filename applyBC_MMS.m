function [L_bc, b_bc] = applyBC_MMS(L, b, dFreedom, elDof)
    L_bc = L;
    b_bc = b;

    nEls = size(dFreedom,1);
    idxR = rightLocalIdx(dFreedom, elDof);

    cL = dFreedom(1,1);              % left endpoint
    cR = dFreedom(nEls, idxR);       % right endpoint (P1/P2-safe)

    [L_bc, b_bc] = pinDirichlet(L_bc, b_bc, cL, 0);
    [L_bc, b_bc] = pinDirichlet(L_bc, b_bc, cR, 0);
end

function [A, b] = pinDirichlet(A, b, dof, val)
    b = b - A(:,dof) * val;
    A(:,dof) = 0;
    A(dof,:) = 0;
    A(dof,dof) = 1;
    b(dof) = val;
end
