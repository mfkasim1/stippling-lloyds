%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to clip 2 polygons with clipping polygon as rectangle.
% The algorithm uses Sutherland-Hodgman algorithm.
% Conditions:
%   * The clipping polygon must be rectangle and the subject polygon can be non-convex.
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
% Input:
%   * xc: array that specifies the x-coordinate of the clipping polygon (px1)
%   * yc: array that specifies the y-coordinate of the clipping polygon (px1)
%   * xs: array that specifies the x-coordinate of the subject polygon (px1)
%   * ys: array that specifies the y-coordinate of the subject polygon (px1)
% Output:
%   * xp: array that specifies the x-coordinate of the resulted polygon (px1)
%   * yp: array that specifies the y-coordinate of the resulted polygon (px1)
% Extended:
%   * can include jacobian vector, so xc, yc, xs, ys, and xp, yp will have dimension of (px3)
%   * the rows are: x : [xA, dxA/dxi, dxA/dyi; xB, dxB/dxi, dxB/dyi; ...]
%   * the rows are: y : [yA, dyA/dxi, dyA/dyi; yB, dyB/dxi, dyB/dyi, ...]

function [xp, yp] = clip_polygons_with_rect(xc, yc, xs, ys)
    if (size(xc,2) == 1)
        ixpos = 1;
        iypos = 2;
    else
        ixpos = 1;
        iypos = 4;
    end
    
    % make the coordinates (x,y) for each polygon
    % if with jacobian, it becomes: [xA, dxA/dxi, dxA/dyi, yA, dyA/dxi, dyA/dyi; ...]
    clippingPolygon = [xc, yc];
    subjectPolygon = [xs, ys];
    
    % get the boundary points
    xcmin = min(xc(:,1)); xcmax = max(xc(:,1));
    ycmin = min(yc(:,1)); ycmax = max(yc(:,1));
    
    outputList = subjectPolygon;
    for (i = [1:size(clippingPolygon,1)])
        % break if there are no point left
        if (length(outputList) == 0) break; end;
        
        % copy the output list and clear it
        inputList = outputList;
        outputList = [];
        
        S = inputList(end,:);
        for (iE = [1:size(inputList,1)])
            E = inputList(iE,:);
            % SEedge = [S; S-E];
            
            % check if E is inside the clipEdge
            if (isInside(xcmin, xcmax, ycmin, ycmax, E, i, ixpos, iypos))
                
                % check if S is not inside the clipEdge
                if (~isInside(xcmin, xcmax, ycmin, ycmax, S, i, ixpos, iypos))
                    
                    % add the intersection from S to E with the clipEdge
                    outputList(end+1,:) = getIntersection(xcmin, xcmax,  ycmin, ycmax, S, S-E, i, ixpos, iypos);
                end
                outputList(end+1,:) = E;
            
            % check if S is inside the clipEdge
            elseif (isInside(xcmin, xcmax, ycmin, ycmax, S, i, ixpos, iypos))
                
                % add the intersection from S to E with the clipEdge
                outputList(end+1,:) = getIntersection(xcmin, xcmax,  ycmin, ycmax, S, S-E, i, ixpos, iypos);
            end
            
            S = E;
        end
    end
    
    if (length(outputList) == 0)
        xp = [];
        yp = [];
    else
        xp = outputList(:,ixpos:iypos-1);
        yp = outputList(:,iypos:end);
        % xp = outputList(:,1)';
        % yp = outputList(:,2)';
    end
end

function ret = isInside(xmin, xmax, ymin, ymax, point, i, ixpos, iypos)
    if (i == 1)
        ret = point(ixpos) >= xmin;
    elseif (i == 2)
        ret = point(iypos) <= ymax;
    elseif (i == 3)
        ret = point(ixpos) <= xmax;
    else
        ret = point(iypos) >= ymin;
    end
end

function intersection = getIntersection(xmin, xmax, ymin, ymax, p0, grad, i, ixpos, iypos)
    if (iypos == 2)
        if (i == 1)
            x = xmin;
            y = p0(iypos) + grad(iypos)/grad(ixpos) * (x - p0(ixpos)); % y0+dy/dx*(x-x0);
        elseif (i == 2)
            y = ymax;
            x = p0(ixpos) + grad(ixpos)/grad(iypos) * (y - p0(iypos)); % x0+dx/dy*(y-y0);
        elseif (i == 3)
            x = xmax;
            y = p0(iypos) + grad(iypos)/grad(ixpos) * (x - p0(ixpos)); % y0+dy/dx*(x-x0);
        else
            y = ymin;
            x = p0(ixpos) + grad(ixpos)/grad(iypos) * (y - p0(iypos)); % x0+dx/dy*(y-y0);
        end
        intersection = [x y];
    else
        if (i == 1)
            x = xmin;
            y = p0(iypos) + grad(iypos)/grad(ixpos) * (x - p0(ixpos)); % y0+dy/dx*(x-x0);
            dxdxi = 0;
            dxdyi = 0;
            dydxi = p0(iypos+1) + (grad(iypos+1)*grad(ixpos) - grad(iypos)*grad(ixpos+1))/grad(ixpos)^2*(x-p0(ixpos)) - grad(iypos)/grad(ixpos)*p0(ixpos+1);
            dydyi = p0(iypos+2) + (grad(iypos+2)*grad(ixpos) - grad(iypos)*grad(ixpos+2))/grad(ixpos)^2*(x-p0(ixpos)) - grad(iypos)/grad(ixpos)*p0(ixpos+2);
        elseif (i == 2)
            y = ymax;
            x = p0(ixpos) + grad(ixpos)/grad(iypos) * (y - p0(iypos)); % x0+dx/dy*(y-y0);
            dydxi = 0;
            dydyi = 0;
            dxdxi = p0(ixpos+1) + (grad(ixpos+1)*grad(iypos) - grad(ixpos)*grad(iypos+1))/grad(iypos)^2*(y-p0(iypos)) - grad(ixpos)/grad(iypos)*p0(iypos+1);
            dxdyi = p0(ixpos+2) + (grad(ixpos+2)*grad(iypos) - grad(ixpos)*grad(iypos+2))/grad(iypos)^2*(y-p0(iypos)) - grad(ixpos)/grad(iypos)*p0(iypos+2);
        elseif (i == 3)
            x = xmax;
            y = p0(iypos) + grad(iypos)/grad(ixpos) * (x - p0(ixpos)); % y0+dy/dx*(x-x0);
            dxdxi = 0;
            dxdyi = 0;
            dydxi = p0(iypos+1) + (grad(iypos+1)*grad(ixpos) - grad(iypos)*grad(ixpos+1))/grad(ixpos)^2*(x-p0(ixpos)) - grad(iypos)/grad(ixpos)*p0(ixpos+1);
            dydyi = p0(iypos+2) + (grad(iypos+2)*grad(ixpos) - grad(iypos)*grad(ixpos+2))/grad(ixpos)^2*(x-p0(ixpos)) - grad(iypos)/grad(ixpos)*p0(ixpos+2);
        else
            y = ymin;
            x = p0(ixpos) + grad(ixpos)/grad(iypos) * (y - p0(iypos)); % x0+dx/dy*(y-y0);
            dydxi = 0;
            dydyi = 0;
            dxdxi = p0(ixpos+1) + (grad(ixpos+1)*grad(iypos) - grad(ixpos)*grad(iypos+1))/grad(iypos)^2*(y-p0(iypos)) - grad(ixpos)/grad(iypos)*p0(iypos+1);
            dxdyi = p0(ixpos+2) + (grad(ixpos+2)*grad(iypos) - grad(ixpos)*grad(iypos+2))/grad(iypos)^2*(y-p0(iypos)) - grad(ixpos)/grad(iypos)*p0(iypos+2);
        end
        intersection = [x dxdxi dxdyi y dydxi dydyi];
    end
end

