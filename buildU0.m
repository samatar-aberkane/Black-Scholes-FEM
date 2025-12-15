function U0 = buildU0(nEls, nodes, connect, elDof, dFreedom, pDeg, pType, u0_fun)
m   = max(dFreedom(:));
U0  = zeros(m,1);
cnt = zeros(m,1);

for e = 1:nEls
    x1 = nodes(connect(e,1));
    x2 = nodes(connect(e,2));

    nLoc = elDof(e);

    xiSamples = linspace(-1, 1, nLoc);

    Nmat = zeros(nLoc, nLoc);
    fvec = zeros(nLoc, 1);

    for k = 1:nLoc
        xi = xiSamples(k);

        [N, ~] = shape(xi, e, pDeg, pType);
        Nmat(k, :) = N(:).';

        Sx = (x1*(1 - xi) + x2*(1 + xi))/2;

        fvec(k) = u0_fun(Sx);
    end

    u_loc = Nmat \ fvec;

    % Scatter-add + count
    for a = 1:nLoc
        gi = dFreedom(e, a);
        U0(gi)  = U0(gi) + u_loc(a);
        cnt(gi) = cnt(gi) + 1;
    end
end

% Average shared DOFs
mask = cnt > 0;
U0(mask) = U0(mask) ./ cnt(mask);
end
