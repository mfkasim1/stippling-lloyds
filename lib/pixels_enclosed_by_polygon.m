%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find pixels which is enclosed by a polygon and crossed by its edges.
% Conditions:
%   * The polygon must be convex
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
%   * Width and height of each pixel must be 1
%   * Bottom left corner is (0,0)
% Input:
%   * x0: the x-coordinates of the polygon's points (nP x 1)
%   * y0: the y-coordinates of the polygon's points (nP x 1)
% Output:
%   * xl: the x-coordinates of the pixels that are crossed by the polygon's edges (numL x 1)
%   * yl: the y-coordinates of the pixels that are crossed by the polygon's edges (numL x 1)
%   * xp: the x-coordinates of the pixels that are fully enclosed by the polygon (numP x 1)
%   * yp: the y-coordinates of the pixels that are fully enclosed by the polygon (numP x 1)

function [xl, yl, xp, yp] = pixels_enclosed_by_polygon(x0, y0)
    % initialise output
    xl = [];
    yl = [];
    xp = [];
    yp = [];
    
    % obtain all pixels crossed by the edges
    for (i = [1:length(x0)])
        % obtain the edge's points
        ip1 = i+1;
        if (i == length(x0)) ip1 = 1; end;
        x1 = x0(i);
        y1 = y0(i);
        x2 = x0(ip1);
        y2 = y0(ip1);
        
        % get the pixels crossed by the edge
        [xline, yline] = pixels_crossed_by_line(x1, y1, x2, y2);
        xl = [xl;xline];
        yl = [yl;yline];
    end
    % remove the duplicates and sort by y
    A = [yl, xl];
    C = unique(A, 'rows');
    yl = C(:,1);
    xl = C(:,2);
    
    % obtain the pixels that are fully enclosed by polygon
    ylmin = min(yl);
    ylmax = max(yl);
    for (y = [ylmin:ylmax])
        xs = xl(yl == y);
        
        % get the index between the edges
        xsdiff = diff(xs);
        idx = find(xsdiff > 1, 1);
        if (length(idx) == 0) continue; end;
        x = [xs(idx)+1:xs(idx)+xsdiff(idx)-1]';
        
        % add the points
        xp = [xp;x];
        yp = [yp;y+zeros(size(x))];
    end
end
