%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find pixels which is crossed by a specific line.
% The difference is that this function use vectorization.
% Input:
%   * x1: the x-coordinate of the initial point 
%   * y1: the y-coordinate of the initial point 
%   * x2: the x-coordinate of the destination point 
%   * y2: the y-coordinate of the destination point 
% Output:
%   * xp: array that specifies the x-coordinate of the pixels (px1)
%   * yp: array that specifies the y-coordinate of the pixels (px1)
% Assumptions:
%   * The line also must have greater displacement in y than x
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1

function [xp, yp] = pixels_crossed_by_specific_line_2(x1, y1, x2, y2)
    
    % initialise output
    xp = [];
    yp = [];
    
    % x2 needs to be larger than x1
    swapped = 0;
    if (x2 < x1)
        xt = x1; yt = y1;
        x1 = x2; y1 = y2;
        x2 = xt; y2 = yt;
        swapped = 1;
    end
    
    % calculate the displacement
    dx = x2 - x1;
    dy = y2 - y1;
    
    % get the sign of the displacement
    sy = (dy >= 0) - (dy < 0);
    
    if (dx == 0)
        % if it's a vertical line, then just calculate xp and yp
        yp = [floor(y1):sy:floor(y2)];
        xp = zeros(size(yp)) + floor(x1);
    else
        % get the coordinates where the line intersect with vertical line
        xintersection = [x1, floor(x1)+1:ceil(x2)-1, x2];
        yintersection = y1 + dy / dx * (xintersection - x1);
        
        % get the pixels of the intersection point
        ixintersection = floor(xintersection);
        iyintersection = floor(yintersection);
        
        % need to repeat ixintersection(1:end-1) by this repetition vector to get xp
        repetition = abs(diff(iyintersection)) + 1;
        
        % repeat ixintersection(1:end-1) by repetition matrix to get xp
        % eg: ixintersection(1:end-1) = [4 5 6]; repetition = [2 3 1];
        % result: [4 4 5 5 5 6];
        cumsumRepetition = cumsum(repetition); % eg: [2 5 6]
        bb = zeros([1 cumsumRepetition(end)]); % eg: [0 0 0 0 0 0]
        bb(cumsumRepetition - repetition + 1) = 1; % eg: bb: [1 0 1 0 0 1]
        cumsumBB = cumsum(bb); % eg: cumsum(bb): [1 1 2 2 2 3]
        xp = ixintersection(cumsumBB); % eg: [4 4 5 5 5 6];
        
        % repeat iyintersection(1:end-1) by repetition matrix with increment 1 to get yp
        % eg: iyintersection(1:end-1) = [4 5 7]; repetition = [2 3 1];
        % result: [4 5, 5 6 7, 7];
        aa = ~bb; % eg: [0 1 0 1 1 0];
        cumsumAA = cumsum(aa); % eg: [0 1 1 2 3 3];
        cc = cumsumAA(aa == 0); % eg: [0 1 3];
        dd = cc(cumsumBB); % eg: [0 0 1 1 1 3];
        deltaDD = cumsumAA - dd; % eg: [0 1 0 1 2 0];
        yp = iyintersection(cumsumBB) + deltaDD*sy;
    end
    
    if (swapped)
        xp = xp(end:-1:1);
        yp = yp(end:-1:1);
    end
    
    xp = xp';
    yp = yp';
end
