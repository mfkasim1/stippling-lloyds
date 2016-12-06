%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find area of convex polygons.
% Conditions:
%   * All polygons must be convex
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
% Input:
%   * x1: (pxN) matrix to specify the x-coordinate of the N polygons
%   * y1: (pxN) matrix to specify the y-coordinate of the N polygons
% Output:
%   * S: (1xN) array that specifies the area of polygons

function S = polyareaconvex(x1, y1)
    if (length(x1) == 0) S = 0; return; end;
    ym1 = [y1(end,:);y1(1:end-1,:)];
    yp1 = [y1(2:end,:);y1(1,:)];
    S = sum(x1 .* ym1 - x1 .* yp1, 1) * .5;
end
