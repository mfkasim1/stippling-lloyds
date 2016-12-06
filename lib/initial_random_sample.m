%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating N random sample with probability distribution valMap using rejection algorithm.
% Input:
%   * N: number of sample
%   * valMap: discrete probability function
% Output:
%   * Px, Py: normalised coordinate of the generated sample (Nx1 each)
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * 1 <= Px < size(valMap,2)+1 and 1 <= Py < size(valMap,1)+1

function [Px,Py] = initial_random_sample(N, valMap)
    [Ny, Nx] = size(valMap);
    
    % normalise the probability distribution function
    valMap = valMap / max(valMap(:));
    
    Px = nan([N,1]);
    Py = nan([N,1]);
    while (any(isnan(Px)))
        nanIdx = find(isnan(Px)); % unfilled position index
        nanN = length(nanIdx); % number of unfilled position
        
        % generate coordinates
        x = rand([nanN,1])*(Nx-1e-3)+1+5e-4;
        y = rand([nanN,1])*(Ny-1e-3)+1+5e-4;
        
        % check the range
        ix = floor(x);
        iy = floor(y);
        inside = (ix >= 1) & (ix < Nx+1) & (iy >= 1) & (iy < Ny+1);
        ix = ix .* inside + (1-inside);
        iy = iy .* inside + (1-inside);
        % disp(any(Nx < ix));
        
        % determine the acceptance of the coordinates
        indices = iy + (ix-1) * Ny;
        prob = valMap(indices);
        r = rand([nanN,1]);
        accept = (r < prob) & inside;
        
        % fill the accepted coordinates
        filledIdx = nanIdx(accept);
        Px(filledIdx) = x(accept);
        Py(filledIdx) = y(accept);
    end
end
