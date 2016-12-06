%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying Lloyd's algorithm for weighted area parallelly.
% Input:
%   * Px0, Py0: list of initial points coordinates (each numPoints x 1)
%   * valMap: value of density per pixel (Ny x Nx)
%   * options:
%       * maxIter: maximum number of iterations (default: 50)
% Output:
%   * Px, Py: list of final points coordinates after applying the Lloyd's algorithm (each numPoints x 1)
%   * Ap: area for each voronoi cell formed by the points (numPoints x 1)

function [Px, Py, Ap] = weighted_lloyds_algorithm(Px0, Py0, valMap, options)
    [Ny, Nx] = size(valMap);
    [X,Y] = meshgrid([1:Nx],[1:Ny]);
    crs = [1, 1; 1, Ny+1; Nx+1, Ny+1; Nx+1, 1];
    
    if (nargin < 4)
        options = struct();
    end
    if ~isfield(options, 'maxIter') options.maxIter = 50; end;
    if ~isfield(options, 'verbose') options.verbose = 0; end;
    
    maxIter = options.maxIter;
    
    Px = Px0;
    Py = Py0;
    for (i = [1:maxIter+1])
        disp(i);
        [V, C] = power_bounded(Px, Py, zeros([size(Px0,1),1]), crs);
        
        Ap = zeros([length(C),1]);
        parfor (j = [1:length(C)])
            % get the centroid
            [xp, yp] = poly2cw(V(C{j}, 1), V(C{j}, 2));
            [A, xc, yc, ~] = poly_pixel_area_cm_inertia(xp, yp, valMap);
            Ap(j) = A;
            
            if (i < maxIter + 1)
                % update the robots positions
                Px(j) = xc;
                Py(j) = yc;
            end
        end
        
        if (options.verbose)
            % close all;
            hold off;
            plot(Px, Py, 'k.', 'markers', 1);
            pause(0.1);
        end
    end
end
