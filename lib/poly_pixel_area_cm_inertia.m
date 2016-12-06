%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the weighted area, centre of mass position, and moment of inertia of a polygon (Px,Py)
% with weight per pixel described in valMap.
% Input:
%   * Px, Py: the normalised coordinate of vertices of the polygon (px1)
%   * valMap: discretised value of weight
%        (i.e. at position [i,i+1) x [j,j+1), the value of weight, f(x), is constant and equal to valMap(i,j))
%        (default: all ones)
% Output:
%   * A: weighted area
%   * xcm, ycm: position of the weighted centre of mass
%   * I0: moment of inertia with respect to the origin (0,0)
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * 1 <= Px < size(valMap,2)+1 and 1 <= Py < size(valMap,1)+1

function [A, xcm, ycm, I0] = poly_pixel_area_cm_inertia(Px, Py, valMap)
    
    % get the pixels enclosed and crossed by the polygon
    [xl, yl, xp, yp] = pixels_enclosed_by_polygon(Px, Py);
    
    A = 0;
    Ax = 0;
    Ay = 0;
    I0 = 0;
    
    % for the pixels enclosed by the polygon
    for (i = [1:length(xp)])
        ix = xp(i);
        iy = yp(i);
        
        % check the range
        if ((iy < 1) || (iy >= size(valMap,1)+1)) continue; end;
        if ((ix < 1) || (ix >= size(valMap,2)+1)) continue; end;
        rho = valMap(iy, ix);
        
        A  = A  + rho;
        Ax = Ax + rho*(ix+.5);
        Ay = Ay + rho*(iy+.5);
        I0 = I0 + rho*(1/6+(ix+.5)^2+(iy+.5)^2);
    end
    
    % for the pixels crossed by the polygon's edge(s)
    for (i = [1:length(xl)])
        ix = xl(i);
        iy = yl(i);
        
        % check the range
        if ((iy < 1) || (iy >= size(valMap,1)+1)) continue; end;
        if ((ix < 1) || (ix >= size(valMap,2)+1)) continue; end;
        rho = valMap(iy, ix);
        
        % construct the voxel coordinate
        xvox = [0;0;1;1] + ix;
        yvox = [0;1;1;0] + iy;
        
        [xpoly, ypoly] = clip_polygons_with_rect(xvox, yvox, Px, Py);
        [dSIntersect, dxIntersect, dyIntersect] = polycmconvex(xpoly, ypoly);
        [Ix, Iy] = polyinertiaconvex(xpoly, ypoly);
        
        A  = A  + rho * dSIntersect;
        Ax = Ax + rho * dSIntersect * dxIntersect;
        Ay = Ay + rho * dSIntersect * dyIntersect;
        I0 = I0 + rho * (Ix + Iy);
    end
    
    if (A ~= 0)
        xcm = Ax / A;
        ycm = Ay / A;
    else
        xcm = mean(Px);
        ycm = mean(Py);
    end
end
