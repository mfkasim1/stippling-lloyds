%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find inertia with respect to origin of convex polygons.
% Conditions:
%   * All polygons must be convex
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
% Input:
%   * x1: (pxN) matrix to specify the x-coordinate of the N polygons
%   * y1: (pxN) matrix to specify the y-coordinate of the N polygons
% Output:
%   * Ix, Iy: moment of inertia with respect to x and y axis

function [Ix, Iy] = polyinertiaconvex(x1, y1)
    if (length(x1) == 0) Ix = 0; Iy = 0; return; end;
    
    xp1 = [x1(2:end,:);x1(1,:)];
    yp1 = [y1(2:end,:);y1(1,:)];
    
    Ix = sum((y1.*y1 + y1.*yp1 + yp1.*yp1) .* (xp1.*y1 - x1.*yp1))/12;
    Iy = sum((x1.*x1 + x1.*xp1 + xp1.*xp1) .* (xp1.*y1 - x1.*yp1))/12;
    
    if (Ix < 0) disp('Ix negative'); end;
end
