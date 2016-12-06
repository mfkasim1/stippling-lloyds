%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find pixels which is crossed by a general line.
% It uses pixels_crossed_by_specific_line algorithm but considering more general cases.
% Input:
%   * x1: the x-coordinate of the initial point 
%   * y1: the y-coordinate of the initial point 
%   * x2: the x-coordinate of the destination point 
%   * y2: the y-coordinate of the destination point 
% Output:
%   * xp: array that specifies the x-coordinate of the pixels (1xp)
%   * yp: array that specifies the y-coordinate of the pixels (1xp)
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1

function [xp, yp] = pixels_crossed_by_line(x1, y1, x2, y2)
    % calculate the displacement
    dy = y2 - y1;
    dx = x2 - x1;
    
    if ((abs(dy) > abs(dx)))
        [xp, yp] = pixels_crossed_by_specific_line_2(x1, y1, x2, y2);
    else
        [yp, xp] = pixels_crossed_by_specific_line_2(y1, x1, y2, x2);
    end
end
