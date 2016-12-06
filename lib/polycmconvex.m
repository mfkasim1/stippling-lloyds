%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find area and centre of mass position of convex polygons.
% Conditions:
%   * All polygons must be convex
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
% Input:
%   * x1: (pxN) matrix to specify the x-coordinate of the N polygons
%   * y1: (pxN) matrix to specify the y-coordinate of the N polygons
% Output:
%   * A: (1xN) array that specifies the area of polygons
%   * xcm, ycm: (1xN) array that specifies the position of centre of mass of polygons

function [A, xcm, ycm] = polycmconvex(x1, y1)
    if (length(x1) == 0) A = 0; xcm = 0; ycm = 0; return; end;
    
    xp1 = [x1(2:end,:);x1(1,:)];
    yp1 = [y1(2:end,:);y1(1,:)];
    
    A = sum(xp1.*y1 - x1.*yp1)*.5;
    
    if (A ~= 0)
        xcm = sum((x1+xp1) .* (xp1.*y1 - x1.*yp1))/(6*A);
        ycm = sum((y1+yp1) .* (xp1.*y1 - x1.*yp1))/(6*A);
    else
        xcm = mean(x1);
        ycm = mean(y1);
    end
    
    if (A < 0) disp('area negative'); end;
    if (xcm < 0) disp('xcm negative'); end;
end
