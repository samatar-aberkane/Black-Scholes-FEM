function plot3D(nEls, nodes, connect, elDof, dFreedom, pDeg, pType, U, time)
% PLOT3D_ELEMENTS
%   Element-colored 3D wireframe of U(tau,x) for 1D FEM.
%
%   U    : m x (Nt+1) matrix of global DOF values (columns = time steps)
%   time : (Nt+1)-vector of times (tau)

    % Sizes
    m   = max(max(dFreedom));      % # global DOFs
    Nt1 = length(time);            % Nt+1

    if size(U,1) ~= m || size(U,2) ~= Nt1
        error('plot3D_elements: U must be m x (Nt+1)');
    end

    % For convenience, local copy of U with time along rows
    % U_t(n,gi) = U(gi,n)
    U_t = U.';                      % size Nt1 x m

    % Plot settings
    nx   = 8;                       % sample points per element for plotting
    cmap = lines(nEls);             % distinct color per element
    lw   = 1.5;

    figure; hold on;

    for ne = 1:nEls
        % Element geometry
        x1 = nodes(connect(ne,1));
        x2 = nodes(connect(ne,2));

        % Local DOF indices
        nLoc = elDof(ne);
        g    = dFreedom(ne,1:nLoc);   % global DOF indices of this element

        % Sample points on this element
        xi = linspace(-1, 1, nx);     % reference coords
        x  = linspace(x1, x2, nx);    % physical coords

        % Precompute shape functions at sample points for this element
        N_all = zeros(nx, nLoc);
        for j = 1:nx
            [N, ~] = shape(xi(j), ne, pDeg, pType);
            N_all(j,:) = N(:).';      % row: 1 x nLoc
        end

        % Color for this element
        col = cmap(ne,:);

        % Time loop: draw polyline for this element at each tau
        for n = 1:Nt1
            % Local DOF values at this time step (row vector 1 x nLoc)
            uLoc = U_t(n, g);         % from global -> local

            % Reconstruct u(x) on sample points
            u_vals = N_all * uLoc.';  % nx x 1

            % Draw 3D line for this element at time tau_n
            t_n = time(n);
            plot3(x, t_n*ones(size(x)), u_vals.', ...
                  'Color', col, 'LineWidth', lw);
        end
    end

    xlabel('S (underlying)');
    ylabel('\tau (time)');
    zlabel('U(\tau,S)');
    title('Element-colored wireframe of FEM solution');

    grid on;
    box on;
    view(135, 30);
    hold off;
end
