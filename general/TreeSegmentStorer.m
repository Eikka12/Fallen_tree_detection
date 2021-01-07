classdef TreeSegmentStorer < handle & matlab.mixin.Copyable
    %TreeSegmentStorer A class for storing and processing tree segments
    
    properties
        trees = {}; 
    end
    
    methods
        function obj = TreeSegmentStorer(treeSegmentPointArray,...
                treeSegmentLineArray)
            %TreeSegmentStorer Construct an instance of this class
            %   obj =
            %   TreeSegmentStorer(treeSegmentPointArray, 
            %   treeSegmentLineArray)
            %   The constructor can be called with either zero or two input
            %   arguments. The optional input arguments are:
            %   treeSegmentPointArray: A cell array in which each cell
            %   contains the points belonging to one tree segment.
            %   treeSegmentLineArray: A nx4 array in which each row
            %   contains the endpoints of one line segment. The format of
            %   the line segment is startX, endX, startY, endY.
            %
            %   The length of treeSegmentPointArray must be the same as the
            %   number of rows in treeSegmentLineArray.
                
            if nargin == 2
                
                if length(treeSegmentPointArray) ~= size(...
                        treeSegmentLineArray,1)
                    error(strcat("The length of the input argument ",...
                        "treeSegmentPointArray must be the same as the",...
                        " number of rows in treeSegmentLineArray"))
                end
                
                % Create TreeSegment objects and store them
                for t = 1:length(treeSegmentPointArray)
                    obj.trees{t} = TreeSegment(treeSegmentPointArray{t},...
                        treeSegmentLineArray(t,:));
                end
                
            elseif nargin ~= 0
                error(strcat("The contructor must be called with 0 or ",...
                    "2 input arguments"))
            end
        end
        
        function addTreeSegments(obj,segments)
            %addTreeSegments Adds the given tree segments to the storer. If
            %there are more than one segments, they must be stored in a row
            %cell array 

            obj.trees = [obj.trees,segments];
        end
        
        function lengths = calculateLengths(obj)
            %calculateLengths Calculates the lengths of the tree segments
            
            lengths = zeros(size(obj.trees));
            for t = 1:length(obj.trees)
                lengths(t) = obj.trees{t}.len;
            end
            
        end
        
        function diameters = calculateDiameters(obj)
            %calculateDiameters Calculates the diameters of the tree
            %segments
            
            diameters = zeros(size(obj.trees));
            for t = 1:length(obj.trees)
                try
                    diameters(t) = obj.trees{t}.diameter;
                % If the tree does not have a diameter
                catch
                    diameters(t) = 0;
                end
            end
        
        end
        
        function lines = getLines(obj)
            %getLines Gets the line segment representations of each tree
            %segment and returns them as a single array
            
            lines = zeros(length(obj.trees),4);
            for t = 1:length(obj.trees)
                lines(t,:) = obj.trees{t}.lineSegment;
            end
        
        end
        
        function segments = getDelineatedSegments(obj)
            %getDelineatedSegments Gets the delineated point cloud
            %representations of each tree segment and returns them in a
            %single cell array. The returned points are in the original
            %coordinate system
            
            segments = cell(length(obj.trees),1);
            for t = 1:length(obj.trees)
                segments{t} = obj.trees{t}.origPts;
            end
        end
        
        function segments = getDelineatedSegmentsTransformed(obj)
            %getDelineatedSegmentsTransformed Gets the delineated point
            %cloud representations of each tree segment and returns them as
            %a single array. The points of each segment are transformed so
            %that the segment lies on the y-axis.
            
            segments = cell(length(obj.trees),1);
            for t = 1:length(obj.trees)
                segments{t} = obj.trees{t}.pts;
            end
        end
        
        function volumes = getVolumes(obj,stepSize,maxDiam)
            %getVolumes Returns the volumes of each tree segment
            %   volumes = getVolumes(obj,stepSize,maxDiam) Takes two input
            %   arguments:
            %   stepSize: The volume of each tree segment is calculated by
            %   dividing the tree segment into intervals with height
            %   stepSize (in meters).
            %   maxDiam: The maximum diameter. If the diameter at interval
            %   HI is over the maximum diameter, this diameter is ignored
            %   and the parts of the tree segment this height interval
            %   divides are merged in the volume calculations. THE
            %   PARAMETER IS OPTIONAL. The default value for maxDiam is 1.
            %
            %   The function returns a vector containing the volumes of
            %   each tree segment (in cubic meters). See
            %   TreeSegment.calculateVolume for more information on how the
            %   volumes are calculated.
            
            if nargin == 2
                maxDiam = 1;
            end
            
            % Create a vector for storing the volumes
            volumes = zeros(size(obj.trees));
            
            % Calculate the volume of each tree segment
            for t = 1:length(obj.trees)
                volumes(t) = obj.trees{t}.calculateVolume(stepSize,maxDiam);
            end
        end
        
        function volume = volumeWithinCircle(obj,circleCenter,radius)
            %volumeWithinCircle Calculates the volume of the sections of
            %all stored tree segments that are located within the given
            %circle.
            %   volume = volumeWithinCircle(obj,circleCenter,radius) Takes
            %   two input arguments:
            %   circleCenter: The coordinates of the centerpoint of a
            %   circle
            %   radius: The radius of the circle
            %
            %   The function calculates the total volume of tree segments
            %   within the given circle. Only the within-circle sections of
            %   the tree segments are included in the volume calculations.
            
            volume = 0;
            % Go through each tree segment
            for t = 1:length(obj.trees)
                % Calculate the volume of the section of the tree segment
                % that is within the sample plot and add it to the total
                % volume
                [~,v] = obj.trees{t}.sectionWithinCircle(circleCenter,...
                    radius);
                volume = volume + v;
            end
        end
        
        % Plotting
        function plotLineSegments(obj,newFig,lineWidth,color)
            %plotLineSegments Plots the stored line segments
            %   plotLineSegments(obj,newFig,lineWidth,color) Takes three
            %   input arguments:
            %   newFig: A logical value determining whether the lines are
            %   plotted on a new figure
            %   lineWidth: The width of the lines
            %   color: A string defining the color of the lines OPTIONAL
            
            % If color argument was not given
            if nargin == 3
                LineProcessor.plotLines(obj.getLines,newFig,lineWidth);
            % If color argument was given
            else
                LineProcessor.plotLines(obj.getLines,newFig,lineWidth,...
                    color);
            end
        end
        
        function plotSegmentsAsPoints(obj,newFig,pointSize,color)
            %plotSegmentsAsPoints Plots the point representations of the
            %stored segments.
            %   plotSegmentsAsPoints(obj,newFig,pointSize,color) Takes
            %   three input arguments:
            %   newFig: A logical value determining whether the lines are
            %   plotted on a new figure
            %   pointSize: The size of the points to be plotted
            %   color: A string defining the color of the points to be
            %   plotted. If this argument is given, all segments are
            %   plotted with the same color. OPTIONAL
            
            if newFig
                figure
                hold on
                axis equal
            end
            
            % Get the point cloud representations of the segments
            segments = obj.getDelineatedSegments;
            
            % If the color argument was not given
            if nargin == 3
                for s = 1:length(segments)
                    pts = segments{s};
                    scatter(pts(:,1),pts(:,2),pointSize,'filled')
                end
                
            % If the color argument was given
            else
                for s = 1:length(segments)
                    pts = segments{s};
                    scatter(pts(:,1),pts(:,2),pointSize,color,'filled')
                end
            end
        end
        
        % Writing segments to a file
        function segmentsToShapefile(obj,filename)
            %segmentsToShapefile Writes the tree segments stored in the
            %object to a single shapefile.
            %   segmentsToShapefile(obj,filename) Takes one input argument:
            %   filename: The name (and path) of the shapefile to be
            %   written.
            %
            %   The function writes the stored tree segments to a
            %   shapefile. The file will consist of polylines each
            %   representing an individual tree segment. Additionally
            %   stored attributes are the length and diameter of the
            %   segment.
            
            % Calculate the lengths and diameters of the segments
            lengths = obj.calculateLengths();
            diameters = obj.calculateDiameters();
            lines = obj.getLines();
            volumes = obj.getVolumes(0.1,0.5); % THESE PARAMETERS CAN BE CHANGED
            
            % Write the segments to a shapefile
            
            % Preallocate
            numSegs = length(obj.trees);
            s = struct('ID',cell(1,numSegs),'Geometry',cell(1,numSegs),...
                'X',cell(1,numSegs),'Y',cell(1,numSegs),'len',...
                cell(1,numSegs),'diameter',cell(1,numSegs),'volume',...
                cell(1,numSegs));
            % Store the information of each segment
            for i = 1:numSegs
                s(i).ID = i;
                s(i).Geometry = 'line';
                s(i).X = lines(i,1:2);
                s(i).Y = lines(i,3:4);
                s(i).len = lengths(i);
                s(i).diameter = diameters(i);
                s(i).volume = volumes(i);
            end

            S = mapshape(s);
            spec = makedbfspec(S);

            shapewrite(S,filename,'DbfSpec',spec)
        end
        
        
            
    end
    
end

