% This script shows a basic example on how to detect fallen trees using the
% function findFallenTrees. Check the documentation of findFallenTrees for
% more information on how to use the function.

%% Add the required filepaths to the MATLAB path
run('../startup.m')

%% Find the fallen trees in the given point cloud
filename = 'ea000001.las';
storer = findFallenTrees(filename);

%% Get the segments
% Get the line segment representations of the detected trees
lines = storer.getLines();
% Get the point cloud representations of the fallen tree segments
segments = storer.getDelineatedSegments();

%% Plot the segments
% Plot the point cloud representations of the segments
newFig = 1; % The segments are plotted on a new figure
pointSize = 5; % The size of plotted points
color = 'r'; % The color of plotted points. OPTIONAL
storer.plotSegmentsAsPoints(newFig,pointSize,color)
% Plot the line representations of the segments
hold on; axis equal;
newFig = 0; % The lines are added to the current figure
lineWidth = 1; % The width of plotted lines
color = 'b'; % The color of plotted lines. OPTIONAL
storer.plotLineSegments(newFig,lineWidth,color)

%% Write the detected fallen trees to a shapefile. 
%The attributes of each segment will contain the estimated length, diameter
%and volume of the segment.
%storer.segmentsToShapefile('test.shp')