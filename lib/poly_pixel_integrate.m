%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the integral \int_P f(x) dx for a given polynomial vertices coordinate (Px, Py)
% with discretise function f(x) equals to valMap.
% Input:
%   * Px, Py: the normalised coordinate of vertices of the polynomial (px1)
%   * valMap: discretised value of f(x)
%        (i.e. at position [i,i+1) x [j,j+1), the value of f(x) is constant and equal to valMap(i,j))
%        (default: all ones)
% Output:
%   * A: integration result
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * 1 <= Px < size(valMap,2)+1 and 1 <= Py < size(valMap,1)+1

function A = poly_pixel_integrate(Px, Py, valMap)
    if (nargin < 3)
        Nx = max(Px);
        Ny = max(Py);
        valMap = ones([Ny Nx]);
    end
    
    % get the pixels enclosed and crossed by the polygon
    [xl, yl, xp, yp] = pixels_enclosed_by_polygon(Px, Py);
    
    A = 0;
    
    % for the pixels enclosed by the polygon
    for (i = [1:length(xp)])
        ix = xp(i);
        iy = yp(i);
        
        % check the range
        if ((iy < 1) || (iy >= size(valMap,1)+1)) continue; end;
        if ((ix < 1) || (ix >= size(valMap,2)+1)) continue; end;
        
        A = A + valMap(iy, ix);
    end
    
    % for the pixels crossed by the polygon's edge(s)
    for (i = [1:length(xl)])
        ix = xl(i);
        iy = yl(i);
        
        % check the range
        if ((iy < 1) || (iy >= size(valMap,1)+1)) continue; end;
        if ((ix < 1) || (ix >= size(valMap,2)+1)) continue; end;
        
        % construct the voxel coordinate
        xvox = [0;0;1;1] + ix;
        yvox = [0;1;1;0] + iy;
        
        [xpoly, ypoly] = clip_polygons_with_rect(xvox, yvox, Px, Py);
        dSIntersect = polyareaconvex(xpoly, ypoly);
        
        A = A + dSIntersect * valMap(iy, ix);
    end
end
