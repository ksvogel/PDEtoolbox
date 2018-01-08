%
% MAKE_MATRIX_RHS_CIRCLEPROBLEM -- make the matrix corresponding to the
%                                  discrete Laplacian on the unit square
%                                  excluding a specified circle inside
%                                  the square f is the rhs corresponding
%                                  to u=1 on the circle equations are
%                                  solved everywhere, inside and out,
%                                  for simplicity
%
%     input:  n   -- number of interior grid points in each direction
%
%     output: A   -- (n*n)x(n*n) matrix of modified discrete Laplacian
%             f   -- (n*n)x1 right-hand-side vector
%             phi -- signed distance function, phi<0 inside circle
%                                              phi>0 outside circle
%
function [A,f,phi]=make_matrix_rhs_circleproblem(n);

    % initialize A as the standard Laplacian (scaled)
    %
    e = ones(n,1);
    L1 = spdiags([e -2*e e],-1:1,n,n);
    I1 = speye(n);
    A  = kron(L1,I1) + kron(I1,L1);

    % initalize f to be zeros
    %
    f = zeros(n*n,1);

    % boundary value assumed to be constant
    %
    Ub=1;

    % form arrays of grid point loctions
    %
    dx = 1/(n+1);
    x = (1:n)*dx;
    [x,y]=ndgrid(x);

    % parameters for the embedded circle
    %
    xc = 0.3;
    yc = 0.4;
    rad = 0.15;

    % compute the signed distance function
    %
    phi = sqrt( (x-xc).^2 + (y-yc).^2 ) - rad;

    % loop over the interior points and flag the irregular points
    %
    IJ = [-1  0;
           1  0;
           0 -1;
           0  1;
         ];

    for j=2:n-1
        for i=2:n-1

            % skip the interior points
            %
            if( phi(i,j) < 0 )
                continue
            end


            % check the neighbors
            %
            for k=1:4

               if( phi(i+IJ(k,1),j+IJ(k,2)) < 0 )

                   % approximate distance to the boundary, scaled by h
                   %
                   alpha = phi(i,j)/( phi(i,j) - phi(i+IJ(k,1),j+IJ(k,2)));

                   % compute the distance to the boundary
                   %
                   kr = sub2ind([n,n],i,j);
                   kc = sub2ind([n,n],i+IJ(k,1),j+IJ(k,2));

                   % adjust RHS
                   %
                   f(kr) = f(kr) - Ub/alpha;

                   % adjust diagional element
                   %
                   A(kr,kr) = A(kr,kr) + 1 - 1/alpha;

                   % adjust off-diagional, and enforce symmetry
                   %
                   A(kr,kc) = 0;
                   A(kc,kr) = 0;

               end

            end

        end
    end
