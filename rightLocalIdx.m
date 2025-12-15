function idxR = rightLocalIdx(dFreedom, elDof)
% Returns the local index in dFreedom(ne,:) that corresponds to the RIGHT endpoint.
% Works for P1 (2 DOFs) and common P2 (3 DOFs) orderings.
%
% P2 common orderings:
%   [L, M, R]  -> right index = 3
%   [L, R, M]  -> right index = 2

nEls = size(dFreedom,1);

% P1 case
if elDof(1) == 2
    idxR = 2;
    return;
end

% P2 needs at least 2 elements to infer ordering
if elDof(1) == 3 && nEls >= 2
    % In a conforming mesh: right endpoint DOF of element 1 equals left endpoint DOF of element 2.
    if dFreedom(1,2) == dFreedom(2,1)
        idxR = 2;   % element DOF order is [L, R, M]
        return;
    elseif dFreedom(1,3) == dFreedom(2,1)
        idxR = 3;   % element DOF order is [L, M, R]
        return;
    end
end

% Fallback: assume last local DOF is right endpoint (may be wrong for some codes)
idxR = elDof(nEls);
end
