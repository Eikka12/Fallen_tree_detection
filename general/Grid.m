classdef Grid < handle & matlab.mixin.Copyable
    %Grid A class representing a raster grid.
    %   The class enables forming a raster from a point cloud. Creating an
    %   instance of this class creates a raster with the defined cell size
    %   and divides the point cloud points into the individual cells of the
    %   grid. The cell values can then be assigned based on various metrics
    %   of the point cloud points. 
    
    properties
        points = []; % The laser points from which the grid is created.
        cellCenters = []; % The center points of the grid cells in real-world coordinates.
        pointsInCells = {}; % A cell array that stores the points grouped into the grid cells.
        values = []; % The value of each grid cell.
        R; % A MapCellsReference object that stores the reference between the grid coordinates (rows,columns) and real-world coordinates.
    end
    
    methods
        % Constructor
        
        function obj = Grid(points,cellSize,parallel)
            %Grid Construct an instance of this class.
            %   obj = Grid(points,cellSize,parallel) Takes three input
            %   arguments:
            %   points: An array in which each row represents one laser
            %   point observation. The first three columns of the matrix
            %   contain the x, y and z coordinates of the laser points,
            %   respectively.
            %   cellSize: The cell size of the grid (in meters).
            %   parallel: A logical value defining whether the grouping of
            %   points into grid cells should be run in parallel using
            %   multiple cores. Use 1 for a point cloud with a large number
            %   of points or a small cell size. The default value is 0.
            %
            %   The grid extent is [min(x),max(x)] in the x direction and
            %   [min(y),max(y] in the y direction.
                
            % if the given points array is empty
            if isempty(points)
                warning('The given points array contains no points.')
                obj.points = points;
                obj.cellCenters = [];
                obj.R = nan;
                obj.pointsInCells = {};
                obj.values = [];
            else
                
                % Handle the parallel argument
                if ~exist('parallel','var')
                    parallel = 0;
                end
                
                % Store the laser point observations
                obj.points = points;

                % Calculate and store the center points of the grid cells.
                % The lower left corner of the grid is at ( min(x),min(y) )
                % and thus the center of the grid cell at the lower left
                % corner is at ( min(x), min(y) ) + cellSize/2. The upper
                % right corner of the grid extends slightly (at most
                % cellSize in each direction) over the maximum extent of
                % the laser points.

                [X,Y] = meshgrid( min(points(:,1)) + cellSize/2 :...
                    cellSize : max(points(:,1)) + cellSize/2 , ...
                    min(points(:,2)) + cellSize/2 : cellSize : ...
                    max(points(:,2)) + cellSize/2);

                % Create a referencing object for storing the reference
                % between the grid coordinates and real-world coordinates
                obj.R = maprefcells([X(1,1) - cellSize/2,...
                    X(1,end) + cellSize/2],...
                    [Y(1,1) - cellSize/2, Y(end,1) + cellSize/2],...
                    [size(Y,1),size(X,2)]);

                % Store the center points of each grid cell in a two-column
                % array
                obj.cellCenters = [X(:),Y(:)];
            
                
                % Groupd the points using the grid cells
                obj.pointsInCells = ...
                    PointCloudProcessor.pointsInCellsTile(...
                    obj.cellCenters,obj.points,obj.getCellSize,...
                    parallel,20);


                % Assign zero as the value of each grid cell
                obj.values = zeros(length(X(:)),1);
            end
        end
        
        % Getters
        
        function numRows = getNumRows(obj)
            %getNumRows Returns the number of rows in the grid.
            %   numRows = getNumRows(obj)
            
            numRows = obj.R.RasterSize(1);
        end
        
        function numCols = getNumCols(obj)
            %getNumCols Returns the number of columns in the grid.
            %   numCols = getNumCols(obj)
            
            numCols = obj.R.RasterSize(2);
        end
        
        function cellSize = getCellSize(obj)
            %getCellSize Returns the cell size of the grid (in meters).
            %   cellSize = getCellSize(obj)
            
            cellSize = obj.R.CellExtentInWorldX;
        end
        
        function [grid,R] = getGrid(obj)
            %createArray Returns the grid as an array.
            %   [grid,R] = getGrid(obj)
            %
            %   Returns the grid as a two-dimensional array in which each
            %   cell is given the value defined by obj.values. In addition,
            %   returns the spatial referencing object of the grid.
            
            if isa(obj.R,'map.rasterref.MapCellsReference') 
                grid = reshape(obj.values,obj.getNumRows,obj.getNumCols);
            
                R = obj.R;
            else
                warning('The grid is empty.')
                grid = [];
                R = nan;
            end
        end
        
        function pts = getRemainingPoints(obj)
            %getRemainingPoints Returns all the laser points in the grid
            %cells as one array.
            %   pts = getRemainingPoints(obj) Takes no input arguments.
            %
            %   The function retuns all the points in obj.pointsInCells in
            %   one array. These points are a subset of obj.points. If no
            %   cell-based filtering has been performed for the points in
            %   each grid cell, the laser points returned by this function
            %   are the same as the points in obj.points. However, the
            %   order of the points in the array may differ.
            
            % Get the cell array containing the points in each grid cell
            ptsCell = obj.pointsInCells;
            
            % Join all points in the cell array to a single array
            pts = obj.pointsToArray(ptsCell);
        end
        
        function pts = getPointsInPolys(obj,polys,allPts)
            %getPointsInPolys Returns the points that are located within
            %the given polygons.
            %   pts = getPointsInPolys(obj,polys,allPts) Takes two input
            %   arguments:
            %   polys: A cell array in which each cell contains a nx2 array
            %   which in turn contains the x and y coordinates of the edge
            %   points of one polygon.
            %   allPts: A logical value defining whether the points are
            %   returned from obj.points or obj.pointsInCells.
            %   obj.pointsInCells will differ if the points have been
            %   filtered somehow. 1 for obj.points and 0 for
            %   obj.pointsInCells. NOTE! IF allPts IS SET TO 0, THE POINTS
            %   WILL BE RETURNED FROM THE CELLS WHOSE CENTERS ARE LOCATED
            %   WITHIN THE GIVEN POLYGON. 
            pts = [];
            if allPts
                % Go through all polygons and find the points that are
                % located within them
                for i = 1:length(polys)
                    % Get the indices of the points that are located within
                    % the given polygon
                    indices = inpolygon(obj.points(:,1),obj.points(:,2),...
                        polys{i}(:,1),polys{i}(:,2));
                    
                    % Append the found points to the pts array
                    pts = [pts;obj.points(indices,:)];
                end
                
            else
                % Go through all polygons and find the cells that are
                % located within them
                for i = 1:length(polys)
                    % Get the indices of the cells that are located within
                    % the given polygon
                    indices = obj.getCellsInPoly(polys{i});
                    
                    % Append the found points to the pts array
                    pts = [pts;...
                        obj.pointsToArray(obj.pointsInCells(indices))];
                end
            end
        end
        
        function indices = getCellsInPoly(obj,poly)
            %getCellsInPoly Returns the indices of the grid cells whose
            %centers are within the given polygon.
            %   indices = getCellsInPoly(obj,poly) Takes one input
            %   argument:
            %   poly: A nx2 matrix containing the x and y coordinates of
            %   the edge points of the polygon.
            %
            %   The function finds the indices of the grid cells whose
            %   center is within or on the border of the given polygon.
            
            indices = inpolygon(obj.cellCenters(:,1),...
                obj.cellCenters(:,2),poly(:,1),poly(:,2));
            
        end
        
        function pts = getPointsInCellLocations(obj,locations)
            %getPointsInCellLocations Returns the points from the given
            %cells.
            %   pts = getPointsInCellLocations(obj,locations) Takes
            %   one input argument:
            %   locations: A nx2 array in which each row contains the x and
            %   y coordinates of one cell location.
            %
            %   The function searches for the points located within the
            %   defined cells and returns them as an array. The cells are
            %   given as the locations of their center points.
            
            % Find the cell indices
            I = ismember(obj.cellCenters,locations,'rows');
            
            % Get the points in each cell
            ptsCell = obj.pointsInCells(I);
            
            % Join all points in the cell array into a single array
            pts = obj.pointsToArray(ptsCell);
        end
        
        % Plotters
            
        function fig = plotGrid(obj,plotType)
            %plotGrid Plots the grid in the given plot type.
            %   fig = plotGrid(obj,plotType) Takes one optional input
            %   argument:
            %   plotType: A string defining the plot type. See the
            %   documentation of mapshow for possible plot types.
            %
            %   The function plots the grid in the given plot type. If no
            %   plot type is given, the grid is plotted as a regular data
            %   grid. The values of each grid cell are acquired from
            %   obj.values. The figure handle is returned as output of the
            %   function.
            
            [grid,r] = obj.getGrid;
            
            fig = figure;
            
            % If no input argument was given
            if nargin == 1
                mapshow(grid,r)
                
            else
                mapshow(grid,r,'DisplayType',plotType)
            end
        end
        
        function fig = plotPoints(obj,all,D2)
            %plotPoints Plots the laser points in the grid.
            %   fig = plotPoints(obj,all,D2) Takes two input arguments:
            %   all: A logical value defining whether all laser points are
            %   plotted or just the ones remaining after possible filtering
            %   operations.
            %   D2: A logical value defining whether the points are plotted
            %   in 2D or 3D. 1 for 2D, 0 for 3D.
            %
            %   The function plots the laser points in obj.points as a
            %   scatter plot. The points are plotted in a new figure. The
            %   figure handle is returned as output of the function.
            
            %Plot all points
            if all
                points = obj.points;
                
            % Plot only the points that remain after filtering    
            else
                points = obj.getRemainingPoints;
            end
            
            fig = figure;
            
            if D2
                scatter(points(:,1),points(:,2),2,'filled')
            else
                scatter3(points(:,1),points(:,2),points(:,3),2,'filled')
            end
            axis equal
            
        end
        
        function fig = plotPointsInCells(obj,cells,D2,newFigure)
            %plotPointsInCells Plots the points located within the given
            %cells.
            %   fig = plotPointsInCells(obj,cells,D2,newFigure) Takes three
            %   input arguments:
            %   cells: A nx2 array containing the center points of the
            %   cells whose points are to be plotted.
            %   D2: A logical value defining whether the points are plotted
            %   in 2D. 1 for 2D, 0 for 3D.
            %   newFigure: A logical value defining whether the points are
            %   plotted on a new figure. 1 for yes, 0 for no.
            %
            %   The function plots the points in the given cells. All
            %   points are plotted with a single colour.
            
            % Get the points
            pts = obj.getPointsInCellLocations(cells);
            
            % If the points should be plotted on a new figure
            if newFigure
                % Create new figure and store its handle
                fig = figure;
            % If the points should be plotted on the current figure
            else
                % Store the handle of the current figure
                fig = gcf;
            end
            
            if D2
                scatter(pts(:,1),pts(:,2),2,'filled')
            else
                scatter3(pts(:,1),pts(:,2),pts(:,3),2,'filled')
            end
            
            axis equal
            
        end
        
        % Cell values
            
        function cellVals = calculateMetric(obj,func)
            %calculateMetric Calculates the metric defined by the given
            %function for each grid cell.
            %   cellVals = calculateMetric(obj,func) Takes one input
            %   argument:
            %   func: An anonymous function
            %
            %   The function uses func to calculate the values of each grid
            %   cell. Func is applied to the laser points of each grid
            %   cell.
            
            cellVals = zeros(length(obj.values),1);
            
            for i = 1:length(cellVals)
                try
                    cellVals(i) = func(obj.pointsInCells{i});
                catch
                    cellVals(i) = 0;
                end
            end
            
            % Replace nans with zero
            cellVals(isnan(cellVals)) = 0;
        end
        
        function means = cellMeanHeight(obj,setAsValue)
            %cellMeanheight Calculates the mean height of the laser points
            %within each cell.
            %   means = cellMeanHeight(obj,setAsValue) Takes one input
            %   argument:
            %   setAsValue: A boolean/numeric value defining whether the
            %   calculated values are set as the values of the grid
            %   (obj.values)
            %
            %   The function calculates the mean height of the laser points
            %   of each grid cell. The means are returned as a vector in
            %   which the value in position i represents the mean value of
            %   the grid cell in position i on obj.cellCenters.
            
            func = @(pts) mean(pts(:,3));
            
            means = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = means;
            end
        end
        
        function vals = cellMeanCol(obj,col,setAsValue)
            %cellMeanCol For each cell, calculates the mean of the values
            %in column col of the laser point observations array.
            %   vals = cellMeanCol(obj,col,setAsValue) Takes two input
            %   arguments:
            %   col: The column in the laser point observations array whose
            %   values are used.
            %   setAsValue: A boolean/numeric value defining whether the
            %   calculates values are set as the values of the grid
            %   (obj.values)
            %
            %   The laser points within each cell are represented as an
            %   array in which the first three columns contain the x-, y-
            %   and z-coordinates of the laser points. The other columns
            %   contain additional information (e.g. the laser point
            %   intensity). The function extracts the values from the colth
            %   column of the laser point array and calculates the mean of
            %   these values for each cell separately.
            
            func = @(pts) mean(pts(:,col));
            
            vals = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = vals;
            end
        end
        
        function counts = cellCount(obj,setAsValue)
            %cellCount Calculates the cell count of each grid cell.
            %   counts = cellCount(obj,setAsValue) Takes one input
            %   argument:
            %   setAsValue: A boolean/numeric value defining whether the
            %   calculated values are set as the values of the grid
            %   (obj.values)
            %
            %   The function calculates the number of laser points falling
            %   within each cell. The cell counts are returned as a vector
            %   in which the value in position i represents the cell count
            %   of the grid cell in position i on obj.cellCenters.
            
            func = @(pts) size(pts,1);

            counts = obj.calculateMetric(func);

            if setAsValue
                obj.values = counts;
            end
        end
        
        function maximums = cellMaxHeight(obj,setAsValue)
            %cellMaxHeight Calculates the maximum laser point height of
            %each grid cell
            %   maximums = cellMaxHeight(obj,setAsValue) Takes one input
            %   argument:
            %   setAsValue: A boolean/numeric value defining whether the
            %   calculated values are set as the values of the grid
            %   (obj.values)
            %
            %   The function calculates the maximum height of the laser
            %   points falling within each cell. The maximums are returned
            %   as a vector in which the value in position i represents the
            %   maximum value of the grid cell in position i on
            %   obj.cellCenters.
            
            func = @(pts) max(pts(:,3));
            
            maximums = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = maximums;
            end 
        end
        
        function minimums = cellMinHeight(obj,setAsValue)
            %cellMinHeight Calculates the minimum laser point height of
            %each grid cell
            %   minimums = cellMinHeight(obj,setAsValue) Takes one input
            %   argument:
            %   setAsValue: A boolean/numeric value edfining whether the
            %   calculated values are set as the values of the grid
            %   (obj.values)
            %
            %   The function calculates the minimum height of the laser
            %   points falling within each cell. The minimums are returned
            %   as a vector in which the value in position i represents the
            %   minimum value of the grid cell in position i on
            %   obj.cellCenters.
            
            func = @(pts) min(pts(:,3));
            
            minimums = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = minimums;
            end 
        end
        
        function medians = cellMedianHeight(obj,setAsValue)
            %cellMedianHeight Calculates the median laser point height of
            %each grid cell.
            %   medians = cellMedianHeight(obj,setAsValue) Takes one input
            %   argument:
            %   setAsValue: A boolean/numeric value defining whether the
            %   calculated values are set as the values of the grid
            %   (obj.values)
            %
            %   The function calculates the median height of the laser
            %   points falling within each cell. The medians are returned
            %   as a vector in which the value in position i represents the
            %   maximum value of the grid cell in position i on
            %   obj.cellCenters.
            
            func = @(pts) median(pts(:,3));
            
            medians = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = medians;
            end
        end
        
        function prctiles = cellPercentile(obj,p,setAsValue)
            %cellPercentile Calculates the pth percentile of heights of the
            %laser points within each grid cell.
            %   prctile = cellPercentile(onj,p,setAsValue) Takes two input
            %   arguments:
            %   p: The percentile (e.g. the value 20 means the 20th
            %   percentile)
            %   setAsValue: A logical value defining whether the calculated
            %   values are set as the values of the grid (obj.values)
            %
            %   The function calculates the pth percentile of the laser
            %   point heights within each cell. The percentiles are
            %   returned as a vector in which the value in position i
            %   represents the percentile of the grid cell in position i on
            %   obj.cellCenters.
            
            func = @(pts) prctile(pts(:,3),p);
            
            prctiles = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = prctiles;
            end
            
        end
        
        function standardDevs = cellHeightSTD(obj,setAsValue)
            %cellHeighSTD Calculates the standard deviations of the laser
            %point heights of each cell.
            %   standardDevs = cellHeightSTD(obj,setAsValue) Takes one
            %   input argument:
            %   setAsValue: A logical value defining whether the calculated
            %   values are set as the values of the grid (obj.values)
            
            func = @(pts) std(pts(:,3));
            
            standardDevs = obj.calculateMetric(func);
            
            if setAsValue
                obj.values = standardDevs;
            end
        end
        
        function correlations = cellCorrelation(obj,cols,setAsValue)
            %cellCorrelation calculates the correlations between the
            %defined columns of the laser points within each cell.
            %   correlations = cellCorrelation(obj,cols,setAsValue) Takes
            %   two input arguments:
            %   cols: A 2x1 vector containing the column numbers for which
            %   the correlation is calculated.
            %   setAsValue: A logical value defining whether the calculated
            %   values are set as the values of the grid (obj.values).
            
            correlations = zeros(length(obj.values),1);
            
            % Go through each cell
            for i = 1:length(obj.pointsInCells)
                pts = obj.pointsInCells{i};
                
                % If the cell contains more than one point
                if size(pts,1) > 1
                    % Calculate correlation matrix
                    corMat = corrcoef(pts(:,cols));
                    % Store the correlation
                    correlations(i) = corMat(2);
                
                % If the cell contains only one point, the correlation is
                % set to one
                else
                    correlations(i) = 1;
                end
            end
            
            if setAsValue
                obj.values = correlations;
            end
        end
        
        % Filtering
        
        function removeConnectedComponents(obj,net,minSize)
            %removeConnectedComponents Removes the connected components in
            %the grid that are classified as the removable class.
            %   removeConnectedComponents(obj,net,minSize) Takes three
            %   input arguments:
            %   net: A trained neural network used for classifying the
            %   connected components.
            %   minSize: No connected components that have an area smaller
            %   than minSize will be removed.
            %
            %   The function assigns the number of points within each grid
            %   cell as the grid cell value. Then, it extracts the
            %   connected components of the grid and calculates various
            %   metrics describing the connected components using
            %   regionprops. Next, the function classifies each connected
            %   component using the input neural network. The connected
            %   components that belong to the class 0 are removed from the
            %   grid, meaning that the laser points falling within these
            %   connected components are removed. However, if the area
            %   (number of cells) of the connected component is smaller
            %   than or equal to minSize, the connected component and its
            %   corresponding points are not removed.
            
            % Assign the number of points within the cell as the cell value
            obj.cellCount(1);
            
            % Get the object as a grid and its corresponding spatial
            % reference object
            [grid,~] = obj.getGrid;
            
            % Calculate the metrics of the connected components
            CC = bwconncomp(grid);
            stats = regionprops(CC,'all');
            
            % Classify each connected component
            classes = round(predictClass(stats,net));
            classes([stats.Area] <= minSize) = 1;
            
            % Set the value of all grid cells that belong to connected
            % components classified as 0 to 0
            removePixelsCellArray = CC.PixelIdxList(~logical(classes));
            removePixels = cat(1,removePixelsCellArray{:});
            grid = ones(size(grid));
            grid(removePixels) = 0;
            
            % Remove the laser points that belong to the connected
            grid = grid(:);
            obj.pointsInCells(grid == 0) = cell(...
                size(obj.pointsInCells(grid == 0)));
        end
        
        % Line Searching
        
        function lines = findLines(obj,angleInterval,numOffsets,...
                minPeakVal,pDist,maxSeparation)
            %findLines Finds lines in the data
            %   lines = findLines(obj,angleInterval,numOffsets,minPeakVal,
            %   pDist,maxSeparation)
            %   Takes 5 input arguments:
            %   angleInterval: The angle bin interval used in the Hough
            %   transformation
            %   numOffsets: The number of offset bins used in the Hough
            %   transformation
            %   minPeakVal: The minimum bin count for the bin to be
            %   detected as a peak
            %   pDist: The maximum distance from a point to the line for
            %   the point to be determined as belonging to the line. Used
            %   when deleting the points after each iteration of the Hough
            %   transformation and when determining the exact line segment
            %   locations.
            %   maxSeparation: The maximum separation between points for
            %   the points to be assigned to the same group. Used when
            %   determining the exact line segment locations.
            %
            %   The function calls the function houghLines to search for
            %   lines in each grid cell. The function houghLines uses the
            %   Hough transformation to find lines in the point cloud. See
            %   the function documentation of houghLines for more detailed
            %   information on how the line search works and what is the
            %   significance of each parameter.
            %
            %   The function returns a cell array. Each cell of the cell
            %   array contains a numeric array that stores the end points
            %   of the lines of a specific grid cell. Each row of the
            %   numeric array represents one line and the columns of the
            %   array represent the start x, end x, start y and end y
            %   coordinates of the line, respectively.
            
            lines = cell(size(obj.pointsInCells));
            for i = 1:length(obj.pointsInCells)
                % If the cell contains points
                if ~isempty(obj.pointsInCells{i})
                    lines{i} = LineSearcher.houghLines(...
                        obj.pointsInCells{i},angleInterval,numOffsets,...
                        minPeakVal,pDist,maxSeparation);
                % If the cell does not contain any points
                else
                    lines{i} = [];
                end
            end
                
        end
        
    end
    
    methods (Access = private)
        
        function pts = pointsToArray(obj,ptsCell)
            %pointsToArray Moves the points separated into different cells
            %of a cell array into a single array.
            %   pts = pointsToArray(ptsCell) Takes one input argument:
            %   ptsCell: A cell array in which each cell contains an array
            %   of laser point observations.
            %
            %   The function joins the arrays of points in the cell array
            %   into a single array of laser points.
            
            % Create an array for storing the points
            pts = zeros(size(obj.points));
            
            % Go through the points in each grid cell
            idx = 1;
            for i = 1:length(ptsCell)
                % Get the number of points in the cell
                n = size(ptsCell{i},1);
                
                % Assign the points in the cell to the array that stores
                % all points
                pts(idx:idx + n-1,:) = ptsCell{i};
                
                % Update the position for assigning the points of the next
                % cell
                idx = idx + n;
                
            end
            
            % Remove the empty rows in pts
            pts = pts(1:idx-1,:);
        end
        
    end
    
end

