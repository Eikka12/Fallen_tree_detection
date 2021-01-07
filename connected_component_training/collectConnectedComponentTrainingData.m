function [DWStats,nonDWStats] = collectConnectedComponentTrainingData(...
    filenames,filter)
%collectConnectedComponentTrainingData A function for manually collecting
%training data for connected component classification.
%   [DWStats,nonDWStats] = collectConnectedComponentTrainingData(...
%   filenames,Name,Value)
%   Creates a manual labelling process for collecting training data of
%   connected components that represent and don't represent fallen trees.
%   Connected components are groups of cells in a binary image in which
%   each cell is in the 8-neighborhood of at least one other cell in the
%   group. The training data can be used for training a connected component
%   classifier that distiguishes connected components belonging to fallen
%   trees from connected components not belonging to fallen trees.
%
%   Input arguments:
%   filenames (REQUIRED): A cell array containing the full or relative
%   paths to the las files from which the connected components will be
%   extracted.
%   classes (NAME-VALUE): A vector containing the classes of laser points
%   to be included in the point cloud. By default, all classes apart from
%   the ground class (2) are included.
%   hRange (NAME-VALUE): A 2-element vector containing the lower and upper
%   limit (in metres) of the height range from which the points are
%   extracted. By default, the height range is [0.2 1].
%
%   Output arguments:
%   DWStats: A structure array containing properties of connected
%   components representing parts of fallen trees.
%   nonDWStats: A structure array containing properties of connected
%   components that do not belong to fallen trees.

arguments
    filenames (1,:) cell
    filter.classes (1,:) double = [1,3:99999];
    filter.hRange (1,2) double = [0.2 1]
end

% Create a PointCloudProcessor object to access required functions
processor = PointCloudProcessor;

% Loop through all given files
for i = 1:length(filenames)
    % Read the file into a structure array
    d = lasdata(filenames{i});
    
    % Extract the given classes, attributes and height range from the
    % point cloud and store the remaining points to a matrix
    M = filterPCloud(d,filter,processor);
    
    % Move on to the next file if the filtered point cloud contains no
    % points
    if isempty(M)
        continue
    end
    
    % Transform the point cloud into a binary image (cell size 0.25). The
    % binary image is a top-view representation of the point cloud. The
    % value of each cell of the image determines whether any points fall
    % within the cell.
    BI = createBinaryImage(M,0.25);
    
    % Find connected components in the binary image and calculate their
    % properties
    stats = findConnectedComponents(BI);
    
    % Ask user to label the connected components of the current file
    [DW,nonDW] = labelComponents(BI);
    
    % Store the properties of the labelled connected components
    % If first file
    if i == 1
        % Create structures for storing the properties of the two different
        % types of connected components. Store the extracted properties
        % from the current file
        DWStats(1:length(DW)) = stats(DW);
        nonDWStats(1:length(nonDW)) = stats(nonDW);
    % If not first file
    else
        % Store the extracted properties from the current file
        DWStats(length(DWStats)+1:length(DWStats)+length(DW)) = stats(DW);
        nonDWStats(length(nonDWStats)+1:length(nonDWStats)+...
            length(nonDW)) = stats(nonDW);
    end
end

end

function M = filterPCloud(d,filter,processor)
% Extracts the given classes and the specified height range from the given
% point cloud. Stores the remaining points to a matrix.

M = processor.pCloud2Matrix(d,'classes',filter.classes);
M = processor.selectHeightRange(M,filter.hRange);

end

function BI = createBinaryImage(M,cellSize)
% Creates a binary image of the given laser points.

raster = Grid(M,cellSize);
raster.cellCount(1);
[BI,~] = raster.getGrid;

end

function stats = findConnectedComponents(BI)
% Finds connected components in the given binary image. Calculates their
% properties.
CC = bwconncomp(BI);
stats = regionprops(CC,'all');

end

function [DW,nonDW] = labelComponents(BI)
% Plots the connected components and asks the user to label the connected
% components as either fallen tree or non fallen tree components.

% Show the binary image as well as a binary image in which each
% connected component has been numbered
figure
subplot(1,2,1)
vislabels(bwlabel(BI))
set(gca,'YDir','normal')
subplot(1,2,2)
imshow(BI)
set(gca,'YDir','normal')

% Ask user to specify the connected components belonging to fallen
% trees and not belonging to fallen trees
DW = str2num(input(strcat("Give the numbers of the connected ",...
    "components belonging to fallen trees separated with a space ",...
    "or comma:\n"),'s'));
nonDW = str2num(input(strcat("Give the numbers of the connected ",...
    "components not belonging to fallen trees separated with a ",...
    "space or comma:\n"),'s'));
close all

end