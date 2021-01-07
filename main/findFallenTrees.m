function storer = findFallenTrees(filename,preprocessing,...
    connCompClassification,detection,delineation,finalClassification)
%findFallenTrees Finds fallen trees from the given point cloud.
%   storer = findFallenTrees(filename,classes,attributes,hRange,
%   connCompClassifier,minSize,detectionCellSize,angleInterval,numOffsets,
%   minPeakVal,pDist,maxSeparation,shortTH,slopeDiff,distDiff,overlapTH,
%   radTH,distTH,finalClassifier,cutoffProb,finalCellSize,
%   useConnCompClassification,useFalseTreeRemoval)
%   
%   The function takes the following input arguments of which 'filename' is
%   required and all other arguments are optional Name-Value arguments:
%   
%   filename: The full filepath to the las file to be processed. The only
%   required argument of the function.
%   classes: A vector containing the laser point classes to be used in line
%   detection. All laser points with a class not included in this vector
%   are removed from the point cloud. By default, only the ground class
%   (class 2) is removed from the point cloud.
%   attributes: A cell array defining the laser point attributes to be
%   stored in addition to the point coordinates. The names of the
%   additional attributes to be included must be provided in lower case
%   letters with spaces replaced by the underscore character _. The laser
%   points will be stored in an array in which each row represents one
%   laser point and the columns contain the attributes of the laser points.
%   The first three columns are reserved for the x-,y- and z-coordinates of
%   the laser points. The possible additional columns contain the
%   additional attributes in the defined order. For example, giving
%   {'intensity','user_data'} as input will include the intensity and user
%   data of each laser point in the fourth and fifth column of the array,
%   respectively. By default, no attributes apart from the coordinates will
%   be included.
%   hRange: A vector containing the lower and upper height bound of the
%   point cloud. During the preprocessing step, all laser points except the
%   laser points falling within these bounds will be removed from the point
%   cloud. The default value is [0.2 1].
%   connCompClassifier: A pattern recognition neural network used for
%   classifying connected components into either representing or not
%   representing fallen trees. The classifier can be created with the
%   function trainConnectedComponentClassifier. If the classifier is not
%   provided as input, the classification is performed using the default
%   classifier.
%   minSize: A minimum size of a connected component (as number of cells)
%   that can be removed. The default value is 3 cells.
%   detectionCellSize: The cell size used for line detection. The point
%   cloud is divided into a grid with the defined cells size. The detection
%   is then performed for the points in each cell separately to reduce the
%   computational cost of the algorithm. The default cell size is 20
%   meters.
%   angleInterval: The angle interval used in iterative Hough transform. If
%   the angle interval is small, points must lie strictly on the same line
%   for the line to be detected, whereas a larger angle interval allows
%   larger deviations from the line. The default value is 1 degree.
%   numOffsets: The number of offset bins used in iterative Hough
%   transform. If the number of offset bins is large, points must lie
%   strictly on the same line for the line to be detected, whereas a small
%   number of offset bins allows larger deviations from the line. The
%   default value is 1000 bins.
%   minPeakVal: The minimum number of points that must lie on the same line
%   for the line to be detected in iterative Hough transform.
%   pDist: The distance threshold used in iterative Hough transform. Points
%   located within pDist from a detected line are determined as belonging
%   to the tree represented by the detected line. Such points are removed
%   to ensure that a different line will be found on the next iteration of
%   Hough transform. The default value is 0.5 meters.
%   maxSeparation: Once the points belonging to a detected line have been
%   extracted using the threshold distance (input argument pDist), the
%   points are projected on the line and grouped. Projected points within a
%   distance defined by maxSeparation from each other are assigned to the
%   same group. The point groups are used for determining the range of the
%   detected line. The result is a line segment whose end points are
%   determined by the range of the largest group of projected points along
%   the line. The default value for maxSeparation is 1 meter.
%   shortTH: The threshold segment length. Once the detected line segments
%   have been merged, the remaining line segments with a length shorter
%   than shortTH are removed from the data. The default value is 2.5 meters
%   slopeDiff: The maximum slope difference (in degrees) for two line
%   segments to be considered belonging to the same line.
%   distDiff: The maximum distance from a start or end point of one line
%   segment to the start or end point of another line segment for the lines
%   to be considered belonging to the same line. The default value is 2
%   meters.
%   overlapTH: The overlap threshold (in percents) used for
%   determining whether two line segments can be merged. Line segments with
%   an overlap percentage above the threshold will not be merged. The
%   default value is 10 %.
%   radTH: One of the two distance thresholds used in delineating the point
%   cloud representation of the tree segments. All points within the
%   distance radTH from a line segment are directly assigned to that
%   segment. The default value is 0.5
%   distTH: One of the two distance thresholds used in delineating the
%   point cloud representation of the tree segments. Points located within
%   distTH from a point already assigned to a tree segment are assigned to
%   the same segment. The default value is 0.2
%   finalClassifier: A convolutional neural network used for classifying
%   delineated tree segments as either true or false tree segments. The
%   classifier can be created with the function trainFinalNetwork. If
%   the classifier is not provided as input, the classification is
%   performed using the default classifier.
%   cutoffProb: The cut-off probability. Segments whose probability of
%   being a tree is predicted to be larger than or equal to the cut-off
%   probability are classified as trees. The smaller the cut-off
%   probability the more segments are classified as trees. The cut-off
%   probability must be in the range [0,1]. The default value is 0.25.
%   finalCellSize: The cell size used for final classification. The cell
%   size must be the same as the cell size used for training the final
%   classifier. The default value is 0.025.
%   useConnCompClassification: A logical value that determines whether the
%   connected component classification step is used. The default value is
%   1.
%   useFalseTreeRemoval: A logical value that determines whether the false
%   tree removal step is used. The default value is 1.
%
%   The function extracts fallen trees from a point cloud. The detection
%   process consists of four steps:
%   1. Reading and preprocessing the data. In this step, the point cloud
%   data are read from a las file and transformed into an array format. In
%   addition, the point cloud is filtered using the input arguments
%   'classes', 'attributes' and 'hRange'.
%   2. Filtering the data using connected component analysis. The function
%   transforms the point cloud into a 2D grid with a cell size of 0.25
%   meters. Then, it calculates various shape descriptors for all connected
%   components in the grid. Finally, the function classifies the connected
%   components as either belonging or not belonging to fallen trees using
%   the provided connected component classifier (input argument
%   'connCompClassifier'). The points belonging to components that do not
%   represent fallen trees are removed from the point cloud. However,
%   points belonging to connected components whose size (in number of grid
%   cells) is smaller than the input argument 'minSize' are not removed.
%   3. Detecting fallen trees as lines in the point cloud. The remaining
%   laser points are divided into a grid with a cell size determined by the
%   input argument 'detectionCellSize'. The function fits lines to the
%   points within each grid cell separately using Hough transform
%   iteratively. The parameters of this line detection process are defined
%   by the input arguments 'angleInterval', 'numOffsets', 'minPeakVal',
%   'pDist' and 'maxSeparation'. The detected line segments are then
%   further processed by merging lines based on the input arguments
%   'slopeDiff', 'distDiff' and 'overlapTH' and removing lines shorter than
%   'shortTH'.
%   4. Delineating fallen tree segments and removing false fallen tree
%   segments. In the last step of the fallen tree detection process, the
%   fallen tree segments are delineated aroung the detected line segments
%   using a region growing procedure that can be tuned with the input
%   arguments 'radTH' and 'distTH'. The delineated line segments are then
%   classified as either true or false fallen tree segments using the
%   provided classifier (input argument 'finalClassifier') and cut-off
%   probability (input argument 'cutoffProb'). The correct cell size used
%   by the classifier must be provided with the input argument
%   'finalCellSize'. Finally, the false fallen tree segments are removed.
%
%   The connected component analysis step (step 2) and false fallen tree
%   removal step (step 4) are optional in the line detection process. The
%   usage of these steps depends on the input arguments
%   'useConnCompClassification' and 'useFalseTreeRemoval'. The function
%   delineates the fallen trees from the point cloud even if the false tree
%   removal step is not used. By default, both steps are used.
%
%   The function returns the detected and delineated fallen tree segments
%   as a TreeSegmentStorer object. The storer contains the line segment
%   representation of each fallen tree segment as well as the points
%   belonging to each fallen tree segment.


arguments
   filename string
   
   preprocessing.classes (1,:) double = []
   preprocessing.attributes (1,:) cell = {}
   preprocessing.hRange (1,2) double...
       {validateHRange(preprocessing.hRange)} = [0.2 1]
   
   connCompClassification.useConnCompClassification...
       {mustBeNumericOrLogical} = 1
   connCompClassification.connCompClassifier = nan
   connCompClassification.minSize (1,1) double {mustBeNonnegative} = 3
   
   detection.detectionCellSize (1,1) double {mustBePositive} = 20
   detection.angleInterval (1,1) double {mustBePositive} = 1
   detection.numOffsets (1,1) double {mustBePositive} = 1000
   detection.minPeakVal (1,1) double {mustBePositive} = 5
   detection.pDist (1,1) double {mustBeNonnegative} = 0.5
   detection.maxSeparation (1,1) double {mustBeNonnegative} = 1
   detection.shortTH (1,1) double {mustBeNonnegative} = 2.5
   detection.slopeDiff (1,1) double {mustBeNonnegative} = 5
   detection.distDiff (1,1) double {mustBeNonnegative} = 2
   detection.overlapTH (1,1) double {mustBeNonnegative} = 10
   
   delineation.radTH (1,1) double {mustBePositive} = 0.5
   delineation.distTH (1,1) double {mustBePositive} = 0.2
   
   finalClassification.useFalseTreeRemoval {mustBeNumericOrLogical} = 1
   finalClassification.finalClassifier = nan
   finalClassification.cutoffProb (1,1) double...
       {mustBeGreaterThanOrEqual(finalClassification.cutoffProb,0),...
       mustBeLessThanOrEqual(finalClassification.cutoffProb,1)} = 0.25
   finalClassification.finalCellSize (1,1) double {mustBePositive} = 0.025
   
end

%%% Handle classifiers
% If no connected component classifier was given and connected component
% classification should be used
if ~isa(connCompClassification.connCompClassifier,'network') &&...
        connCompClassification.useConnCompClassification
    % Use default connected component classifier
    cClassifier = load('connCompClassifier.mat');
    connCompClassification.connCompClassifier =...
        cClassifier.connCompClassifier;
end

% If no final classifier was given and final classification should be used
if ~isa(finalClassification.finalClassifier,'SeriesNetwork') &&...
        finalClassification.useFalseTreeRemoval
    % Use default final classifier
    fClassifier = load('finalClassifier.mat');
    finalClassification.finalClassifier = fClassifier.finalClassifier;
    finalClassification.finalCellSize = 0.025;
end

%% 1. Preprocessing
processor = PointCloudProcessor;

% Read the point cloud data
d = lasdata(filename);

%%% Store the laser points into an array
% If no classes were given, use the default selection (only ground class is
% removed)
if isempty(preprocessing.classes)
    classes = unique(d.get_classification);
    preprocessing.classes = classes(classes ~= 2);
end

M = processor.pCloud2Matrix(d,'classes',preprocessing.classes,...
    'attributes',preprocessing.attributes);

%%% Extract laser points falling within the specified height range
M = processor.selectHeightRange(M,preprocessing.hRange);

% If the point cloud contains no points
if isempty(M)
    % Return an TreeSegmentStorer containing no segments
    storer = TreeSegmentStorer({},[]);
    return
end

%% 2. Connected component filtering

% If the connected component classification should be used
if connCompClassification.useConnCompClassification
    % Create a grid from the point cloud
    raster = Grid(M,0.25);
    % Remove connected components that do not represent fallen trees
    raster.removeConnectedComponents(...
        connCompClassification.connCompClassifier,...
        connCompClassification.minSize)
    M = raster.getRemainingPoints;
end

%% 3. Line detection

% Create a grid from the point cloud
raster = Grid(M,detection.detectionCellSize);
% Find line segments in the point cloud
lines = raster.findLines(detection.angleInterval,detection.numOffsets,...
    detection.minPeakVal,detection.pDist,detection.maxSeparation);
% Merge and remove line segments
lines = LineSearcher.lineCleaning(lines,detection.shortTH,...
    detection.slopeDiff,detection.distDiff,detection.overlapTH);

% If no lines were found
if isempty(lines)
    % Return an TreeSegmentStorer containing no segments
    storer = TreeSegmentStorer({},[]);
    return
end

%% 4. Delineation and final classification
delineator = TreeDelineator(lines,M);
% Delineate the trees
delineator.delineateTrees(delineation.radTH,delineation.distTH);

% If the false tree removal should be used
if finalClassification.useFalseTreeRemoval
    % Classify the delineated trees and remove false trees
    delineator.removeFalseTrees(finalClassification.finalClassifier,...
        finalClassification.finalCellSize,finalClassification.cutoffProb)
end

%% Store the detected and delineated tree segments
storer = TreeSegmentStorer(delineator.trees,delineator.lines);

end

function validateHRange(hRange)

if diff(hRange) < 0
    eidType = 'firstElementMustBeSmaller:firstElementNotSmaller';
    msgType = 'The lower bound must be smaller than the upper bound.';
    throwAsCaller(MException(eidType,msgType))
end

end

