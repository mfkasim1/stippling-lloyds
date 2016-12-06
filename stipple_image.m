%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stippling an image using weighted Lloyd's algorithm.
% The function receives input of the image and it produces 2 svg files of the stipple images: one in black and white and one in colors.
% Inputs:
%   * fdir: directory of the image you want to stipple (it also serves as directory where the results are saved)
%   * fname: name of the image file you want to stipple
%   * N: number of dots in the stipple image (optional, default: 20000)
% Outcomes:
%   * drawings.svg: stippled image file in color
%   * drawings-bw.svg: stippled image file in black and white
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stipple_image(fdir, fname, N)
    addpath('lib');
    if (nargin < 3)
        N = 20000;
    end
    
    % normalising the directory
    fdir = strrep(fdir, '\', '/');
    if (fdir(end) ~= '/')
        fdir = strcat(fdir, '/');
    end
    
    % process the image
    colorImg = imread(strcat(fdir, fname));
    img = double(rgb2gray(colorImg));
    img = 1 - img / max(img(:));

    % deploy initial random samples using a simple rejection method
    [Px0, Py0] = initial_random_sample(N, img);

    % use lloyds' algorithm
    options = struct();
    options.maxIter = 50;
    options.verbose = 1;
    [Px, Py, Ap] = weighted_lloyds_algorithm(Px0, Py0, img, options);
    % save(strcat(fdir, 'results.mat'), 'Px', 'Py', 'Ap', 'img');

    % give color to the points
    disp('Obtaining the color of each point');
    [Ny, Nx] = size(img);
    crs = [1, 1; 1, Ny+1; Nx+1, Ny+1; Nx+1, 1];
    [V, C] = power_bounded(Px, Py, zeros([size(Px0,1),1]), crs);
    colorImg = double(colorImg);
    whites = ones(size(colorImg));
    % rgbs = zeros(size(Px,1), 3);
    reds = zeros(size(Px));
    greens = zeros(size(Px));
    blues = zeros(size(Px));
    parfor (j = [1:size(Px,1)])
        disp(j);
        [xp, yp] = poly2cw(V(C{j}, 1), V(C{j}, 2));
        [Ared, ~, ~, ~] = poly_pixel_area_cm_inertia(xp, yp, colorImg(:,:,1));
        [Agreen, ~, ~, ~] = poly_pixel_area_cm_inertia(xp, yp, colorImg(:,:,2));
        [Ablue, ~, ~, ~] = poly_pixel_area_cm_inertia(xp, yp, colorImg(:,:,3));
        [A, ~, ~, ~] = poly_pixel_area_cm_inertia(xp, yp, whites);
        reds(j) = Ared / A;
        greens(j) = Agreen / A;
        blues(j) = Ablue / A;
    end
    rgbs = round([reds, greens, blues]);
    % save(strcat(fdir, 'colors.mat'), 'rgbs');

    % write a csv file of the dots
    disp('Writing the svg file');
    svgHeight = Ny+1;
    svgWidth = Nx+1;
    dotRadius = sqrt(Nx*Ny / N) / sqrt(max(img(:)) / mean(img(:))) / 2.5;
    svgOpening = sprintf('<svg height="%d" width="%d">', svgHeight, svgWidth);
    svgClosing = '</svg>';
    svgDots = '';
    svgDotsBW = '';
    parfor (j = [1:size(Px,1)])
        disp(j);
        svgDot = sprintf('<circle cx="%f" cy="%f" r="%f" stroke-width="0" fill="#%02x%02x%02x" />', Px(j), Py(j), dotRadius, rgbs(j,1), rgbs(j,2), rgbs(j,3));
        svgDotBW = sprintf('<circle cx="%f" cy="%f" r="%f" stroke-width="0" fill="#000000" />', Px(j), Py(j), dotRadius);
        svgDots = strcat(svgDots, svgDot);
        svgDotsBW = strcat(svgDotsBW, svgDotBW);
    end
    svgStr = strcat(svgOpening, svgDots, svgClosing);
    svgStrBW = strcat(svgOpening, svgDotsBW, svgClosing);

    fileID = fopen(strcat(fdir, 'drawings.svg'),'w');
    fprintf(fileID, svgStr);
    fclose(fileID);

    fileID = fopen(strcat(fdir, 'drawings-bw.svg'),'w');
    fprintf(fileID, svgStrBW);
    fclose(fileID);
end
