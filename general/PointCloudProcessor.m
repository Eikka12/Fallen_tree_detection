classdef PointCloudProcessor
    %PointCloudProcessor A class that contains methods for processing laser
    %point clouds.
    
    methods
        function obj = PointCloudProcessor()
            %PointCloudProcessor Construct an instance of class
            %PointCloudProcessor
            %   obj = PointCloudProcessor() Takes no input arguments.
        end

        function M = pCloud2Matrix(obj,pCloud,varargin)
            %pCloud2Matrix Transforms the data in a lasdata object into
            %matrix format.
            %   M = pCloud2Matrix(pCloud,varargin) Takes the following
            %   input arguments:
            %   pCloud = A lasdata object REQUIRED
            %   classes = A vector containing the classes whose points will
            %   be stored. OPTIONAL NAME-VALUE
            %   attributes = A cell array containing the additional
            %   attributes that will be stored. When defining the
            %   attributes, spaces should be replaced with a _ character.
            %   For example the user data attribute should be defined as
            %   user_data (see the getter methods of lasdata objects for
            %   the possible attributes that can be stored. OPTIONAL
            %   NAME-VALUE
            %
            %   The function converts the lasdata object into a matrix in
            %   which each row represents one laser point observation. The
            %   function can extract only points belonging to a certain
            %   class (the classes to be extracted are defined in the
            %   classes input argument) or all laser points in the lasdata
            %   object. The x, y and z (h) coordinate of the laser points
            %   will always be extracted and will be stored in the first,
            %   second and third columns of the returned matrix,
            %   respectively. The function can extract additional
            %   attributes, which should be defined in the attributes input
            %   argument. The attributes will be stored in order starting
            %   from the fourth column of the returned matrix.
            %
            %   Examples:
            %   pCloud2Matrix(pCloud) returns a matrix containing the x, y
            %   and z (h) coordinates of all observations in the lasdata
            %   object (pCloud).
            %   pCloud2Matrix(pCloud,'classes',[1 2],
            %   'attributes',{'intensity'}) returns a matrix containing
            %   the x, y, and z (h) coordinates as well as the intensity of
            %   the laser points belonging to classes 1 and 2.
            
            % Parse inputs
            p = inputParser;
            
            checkPCloud = @(x) isequal(class(x),'lasdata');
            addRequired(p,'pCloud',checkPCloud)
            addParameter(p,'classes',@isnumeric)
            addParameter(p,'attributes',@iscell)
            
            parse(p,pCloud,varargin{:})
            
            % Create a matrix containing the x, y and z coordinates of all
            % laser points
            M = [pCloud.x pCloud.y pCloud.z];
            
            % Process the possible optional attributes
            if ~isempty(varargin)
                if length(varargin) == 4
                    if isequal(varargin{1},'attributes')
                        M = obj.addAttributes(M,pCloud,varargin{2});
                        M = obj.getPointsWithClasses(M,pCloud,varargin{4});
                    else 
                        M = obj.addAttributes(M,pCloud,varargin{4});
                        M = obj.getPointsWithClasses(M,pCloud,varargin{2});
                    end
                else
                    if isequal(varargin{1},'attributes')
                        M = obj.addAttributes(M,pCloud,varargin{2});
                    else
                        M = obj.getPointsWithClasses(M,pCloud,varargin{2});
                    end
                end
            end
            
        end
        
        function M = selectHeightRange(obj,M,HRange)
            %selectHeightRange Selects the laser points that fall within
            %the given height range.
            %   M = selectHeightRange(M,HRange) Takes two input
            %   arguments:
            %   M: A matrix containing laser point observations.
            %   HRange: [minH,maxH] - The height range from which the laser
            %   points will be extracted.
            %
            %   The function returns the laser points from the given matrix
            %   that fall within the given height range. The height range
            %   is exclusive from the lower bound and inclusive from the
            %   upper bound (]minH,maxH]).
            
            func = @(pts) pts(:,3) > HRange(1) & M(:,3) <= HRange(2);
            
            M = obj.selectPoints(M,func);
           
        end
        
    end
    
    methods (Static)
        
        function pointsCell = pointsInCellsTile(cells,pts,cellSize,...
                parallel,tileSize)
            %pointsInCellsTile Gets the points falling within each cell.
            %The processing is performed in tiles.
            %   pointsCell = pointsInCellsTile(cells,pts,cellSize,parallel,
            %   tileSize)
            %   Takes five input arguments:
            %   cells: A nx2 array in which each row represents the x and y
            %   coordinates of the center point of one cell.
            %   pts: An array containing laser point observations. The
            %   first two columns of the array represent the x and y
            %   coordinates of the laser point observations.
            %   cellSize: The cell size in meters.
            %   parallel: A logical value that determines whether the
            %   points in each cell are searched in parallel or not. Use
            %   parallel search for large datasets or small cell sizes. 0
            %   or false means that the search is not run in parallel.
            %   tileSize: The tile size in meters. The grid and laser
            %   points are separated into tiles whose tile size in the
            %   direction of the x-axis is tilesize. No tiling is done in
            %   the y-direction. The search is performed to each tile
            %   separately to improve the computing time.
            %
            %   The function searches for all laser points falling within
            %   each cell and returns them as a cell array. The laser
            %   points in position i in the cell array belong to the grid
            %   cell in position i in the input argument cells.
            
            % If the cell size is larger than the tile size
            if cellSize > tileSize
                warning(strcat("The cell size is larger than the tile ",...
                "size. The tile size is set equal to the cell size"))
                tileSize = cellSize;
            end
            
            % Get the unique x-coordinates
            cols = unique(cells(:,1));
            
            % Modify the tile size so that it is an equal multiple of the
            % cell size
            modulo = mod(tileSize,cellSize);
            tileSize = tileSize - modulo;
            
            % Define the border locations of each tile
            tileBorders = min(cols) - cellSize/2 : tileSize : max(cols)...
                + cellSize/2;
            
            % Store the data in each tile separately
            
            cellsTiled = cell(1,length(tileBorders));
            ptsTiled = cell(1,length(tileBorders));
            
            % Go through each tile
            for i = 2:length(tileBorders)
                % Cells in tile
                tileLocs = cells(:,1) >= tileBorders(i-1) &...
                    cells(:,1) < tileBorders(i);
                cellsTiled{i-1} = cells(tileLocs,:);
                % Points in tile
                ptsLocs = pts(:,1) >= tileBorders(i-1) &...
                    pts(:,1) < tileBorders(i);
                ptsTiled{i-1} = pts(ptsLocs,:);
            end
            % Store the cells and points that have a larger x-coordinate
            % than the last border
            cellsTiled{end} = cells(cells(:,1) >= max(tileBorders),:);
            ptsTiled{end} = pts(pts(:,1) >= max(tileBorders),:);
            
            % Find the points in each cell
            processor = PointCloudProcessor;
            pointsCell = {};
            % Go through all tiles
            for i = 1:length(cellsTiled)
                pointsCell = [pointsCell,processor.pointsInCells(...
                    cellsTiled{i},ptsTiled{i},cellSize,parallel)];
            end
        end    
        
        function pointsCell = pointsInCells(cells,pts,cellSize,parallel)
            %pointsInCells Gets the points falling within each cell
            %   pointsCell = pointsInCells(cells,pts,cellSize,parallel)
            %   Takes four input arguments:
            %   cells: A nx2 array in which each row represents the x and y
            %   coordinates of the center point of one cell.
            %   pts: An array containing laser point observations. The
            %   first two columns of the array represent the x and y
            %   coordinates of the laser point observations.
            %   cellSize: The cell size in meters.
            %   parallel: A logical value that determines whether the
            %   points in each cell are searched using parallel processing
            %   or not. Use parallel search for large datasets or very
            %   small cell sizes. 0 or false means that the search is not
            %   run in parallel.
            %
            %   The function searches for all laser points falling within
            %   each cell and returns them as a cell array. The laser
            %   points in position i in the cell array belong to the grid
            %   cell in position i in the input argument cells.

            % Cell array for storing the points
            pointsCell = cell(1,size(cells,1));
            % Store the x and y coordinates of the points
            pts2D = pts(:,1:2);
            % The upper and right boundaries of the cells
            cellsmax = cells+cellSize/2;
            % The lower and left boundaries of the cells
            cellsmin = cells-cellSize/2;

            if parallel
                parfor i=1:size(cells,1)
                    % Get a logical vector containing the locations of the
                    % laser points within the cell
                    logVec = PointCloudProcessor.pointsInCell(...
                        cellsmin(i,:),cellsmax(i,:),pts2D);
                    % Get the laser points within the cell and store them
                    pointsCell{i} = pts(logVec,:);
                end
                
            else
                for i=1:size(cells,1)
                    % Get a logical vector containing the locations of the
                    % laser points within the cell
                    logVec = PointCloudProcessor.pointsInCell(...
                        cellsmin(i,:),cellsmax(i,:),pts2D);
                    % Get the laser points within the cell and store them
                    pointsCell{i} = pts(logVec,:);
                end
            end

        end
        
        function logVec = pointsInCell(cellMin,cellMax,pts2D)          
            %pointsInCell Returns a logical vector that determines the
            %positions of the laser points located within the given cell.
            %   logVec = pointsInCell(cellMin,cellMax,pts2D) Takes three
            %   input arguments:
            %   cellmin: A nx2 array containing the lower and left
            %   boundaries of the cells.
            %   cellMax: A nx2 array containing the upper and right
            %   boundaries of the cells.
            %   pts2D: An array containing the x and y coordinates of laser
            %   points.
            %
            %   The function searches for all laser points in pts2D that
            %   are located within the cell defined by cellMin and cellMax.
            %   The positions of these points in pts2D are returned as a
            %   logical vector. Points located on the lower or left
            %   boundary of the cell are included, whereas the points
            %   located on the upper or right boundary of the cell are
            %   excluded.
            
            logVec = pts2D(:,1) < cellMax(1) &...
                pts2D(:,1) >= cellMin(1) & pts2D(:,2) < cellMax(2) &...
                pts2D(:,2) >= cellMin(2);
        end
        
        function pts = pointsInPolygons(points,polygons,separate)
            %pointsInPolygons Returns the laser points falling within the
            %given polygons.
            %   pts = pointsInPolygons(points,polygons,separate) Takes
            %   three input arguments:
            %   points: A matrix containing laser point observations. The
            %   first two columns contain the x and y coordinates of the
            %   observations.
            %   polygons: A cell matrix containing the edge coordinates of
            %   the polygons.
            %   separate: A logical value defining whether the points found
            %   are returned separated by polygon or all together. A
            %   logical true returns a cell array in which each cell
            %   contains the points within one polygon. A logical false
            %   returns all points as a matrix.
            %   
            %   The function returns the laser points that fall within or
            %   on the border of the given polygons.
            
            % Store the points within each polygon separately
            if separate
                
                pts = cell(1,length(polygons));
                for p = 1:length(polygons)
                    pts{p} = PointCloudProcessor.pointsInPolygon(...
                        points,polygons{p});
                end
                    
            % Store the points within each polygon together    
            else
                
                pts = [];
                for p = 1:length(polygons)
                    pts = [pts ; PointCloudProcessor.pointsInPolygon(...
                        points,polygons{p})];
                end
            end
        end
        
        function pts = pointsInPolygon(points,polygon)
            %pointsInPolygon Returns the laser points falling within the
            %given polygon.
            %   pts = pointsInPolygon(points,polygon) Takes two input
            %   arguments:
            %   points: A matrix containing laser point observations. The
            %   first two columns contain the x and y coordinates of
            %   the observations.
            %   polygon: A nx2 matrix containing the edge coordinates of
            %   the polygon.
            %
            %   The function returns the laser points that fall within or
            %   on the border of the given polygon.
            
            in = inpolygon(points(:,1),points(:,2),polygon(:,1),...
                polygon(:,2));
            
            pts = points(in,:);
        end
        
        function M = selectPoints(M,condFunc)
            %selectPoints Selects the points that fulfill the given
            %condition.
            %   M = selectPoints(M,condFunc) Takes two input arguments:
            %   M: An array of laser point observations
            %   condFunc: A function defining the condition used for
            %   selecting the points.
            %
            %   The function returns the points of M that fulfill the given
            %   condition.
            
            % Get the indices of the laser points that fulfill the
            % condition
            indices = condFunc(M);
            
            M = M(indices,:);
            
        end
        
        function dists = distToPoint2D(point,points)
            %distToPoint2D Calculates the 2D Euclidean distance from
            %several points to one reference point.
            %   dists = distToPoint2D(point,points) Takes two input
            %   arguments:
            %   point: A 2-dimensional vector containing the x and y
            %   coordinates of the reference point.
            %   points: A nx2 array containing the x and y coordinates of
            %   the points whose distance to the reference point will be
            %   calculated.
            %
            %   The function calculates the 2D Euclidean distance from each
            %   point in points to the reference point and returns the
            %   distances as a nx1 vector.
            
            dists = sqrt( (points(:,1) - point(1)).^2 +...
                (points(:,2) - point(2)).^2 ); 
        end
        
        function pDensity = pointDensity(M)
            % pointDensity Calculates the point density of the given laser
            % points.
            %   pDensity = pointDensity(M)Takes one input argument:
            %   M: An array of laser points
            %
            %   The function calculates the area of the laser points as the
            %   horizontal convex hull of the points. It then calculates
            %   the point density as the number of laser points divided by
            %   the convex hull area
            
            % Calculate the 2D convex hull of the input points
            k = convhull(M(:,1:2));
            
            % Calculate the area of the convex hull
            area = polyshape(M(k,1:2)).area;
            
            % Calculate the point density
            pDensity = size(M,1)/area;
     
        end
        
    end
       
    methods (Static, Access = private)
        
        function M = getPointsWithClasses(M,pCloud,classes)   
            %getPointsWithClasses Returns the points of input matrix m that
            %belong to the given classes.
            %   M = getPointsWithClasses(M,pCloud,classes) Takes
            %   three input arguments:
            %   M = A matrix containing laser point observations.
            %   pCloud = The lasdata object from which the laser points in
            %   M have been extracted.
            %   classes = A vector containing the classes whose point will
            %   be extracted.
            
            % Get the classes of each laser observation 
            classification = pCloud.get_classification;
            
            % Get the rows that belong to the given classes
            validRows = ismember(classification,classes);
            M = M(validRows,:);
            
        end
        
        function M = addAttributes(M,pCloud,attributes)
            %addAttributes Adds the defined additional attributes to the
            %point cloud matrix.
            %   M = addAttributes(M,pCloud,attributes) Takes three input
            %   arguments:
            %   M = A matrix containing the x, y and z (h) coordinates of
            %   laser point observations.
            %   pCloud = The lasdata object from which the laser points in
            %   M have been extracted.
            %   attributes = A cell array containing the names of the
            %   attributes to be extracted.
            for i = 1:length(attributes)
                % Try to get the attribute data
                try
                    atbs = eval(strcat('pCloud.get_',attributes{i}));
                % If the attribute is not valid
                catch
                    error(strcat("The attribute ",attributes{i},...
                        " is not valid."))
                end
                
                % Check whether attribute data was available
                if length(atbs) == size(M,1)
                    % Add the attribute data to the matrix
                    M = [M double(atbs)];
                else
                    warning(strcat("The attribute ",attributes{i},...
                        " was not added to the matrix,",...
                        " as no data of this attribute was available."))
                end
            end
            
        end
        
    end
    
end

