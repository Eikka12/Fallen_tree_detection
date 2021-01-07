function [segments,labels] = collectTrainingSegments(filenames,cellSize)
%collectTrainingSegments A function for collecting training segments for
%fallen tree segment classification.
%   [segments,labels] = collectTrainingSegments(filenames,cellSize) Takes
%   two input arguments:
%   filenames: A cell array containing full filepaths of the las files from
%   which the training segments are collected.
%   cellSize: The cell size of the binary image that is shown to the user 
%
%   The function detects and delineates fallen tree segments from the input
%   point clouds using iterative Hough transform and region growing. These
%   segments are then shown (as a point cloud and a binary image) to the
%   user and the user is asked to label each segment manually.
%
%   The function returns a cell array containing the labelled segments (as
%   arrays of laser points) and a logical vector containing the label of
%   each segment.

% Create a cell array for storing the training segments and a vector for
% storing the labels
segments = {};
labels = [];

% Go through all provided las files
for f = 1:length(filenames)
    %% Extract tree segments from the point cloud
    % Delineate the tree segments
    storer = findFallenTrees(filenames{f},'useFalseTreeRemoval',0);
   
    % Get the tree segments with points transformed to lie on the y-axis
    segs = storer.getDelineatedSegmentsTransformed();
    
    % If some segments were found in the current point cloud
    if ~isempty(segs)
        
        %% Label and store each tree segment
        for s = 1:length(segs)
            seg = segs{s};

            % If the current segment contains points
            if ~isempty(seg)
                % Store the current segment
                segments{end+1} = seg;
                % Label the current segment
                labels(end+1) = labelSegment(seg,cellSize);
            end
        end
    end
end

% Create a logical vector from the labels
labels = logical(labels);

end

function label = labelSegment(seg,cellSize)
%labelSegment Plots the given segment and asks the user to label it.

% Create a binary image of the segment
raster = Grid(seg,cellSize);
raster.cellCount(1);
[g,~] = raster.getGrid;

% Show the segment as a point cloud and binary image
subplot(1,2,1)
scatter(seg(:,1),seg(:,2),10,'filled')
axis equal
subplot(1,2,2)
imshow(flipud(g))
set(gcf,'Position',get(0,'ScreenSize'));

% Ask user if current segment represents a fallen tree
flag = 1;
while flag 
    ip = input('Does the image represent a tree (1 = yes, 0 = no)? ');
    % If valid input was given
    if ~isempty(ip) && (ip == 0 || ip == 1)
        label = ip;
        flag = 0;
    else
        disp("Input is not valid.")
    end
end

close all;

end

