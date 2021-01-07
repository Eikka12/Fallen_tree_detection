classdef TreeDelineator < handle
    %TreeDelineator Delineates fallen trees around found lines in point
    %cloud
    
    properties
        lines; % Lines
        slopes; % The line slopes
        intercepts; % The line intercepts
        M; % Laser points
        trees; % A cell array in which each cell contains the points of one tree segment
    end
    
    methods
        function obj = TreeDelineator(lines,M)
            %treeDelineator Construct an instance of this class
            %   obj = TreeDelineator(lines,M) Takes two input arguments:
            %   lines: An array of detected lines. Each row of the array
            %   represents the end points of one line. The columns of the
            %   array are startX, endX, startY, endY.
            %   M: An array of laser points. Each row of the array
            %   represents one laser point. The first three columns of the
            %   array contain the x, y and z coordinates of the points,
            %   respectively.
            
            % Order the lines in descending order based on line length
            lens = LineProcessor.distPts2Pts(lines(:,[1 3]),...
                lines(:,[2,4]));
            lines = sortrows([lens,lines],'descend');
            lines(:,1) = [];
            obj.lines = lines;
            
            % Calculate the line slopes and intercepts
            [obj.slopes,obj.intercepts] =...
                LineProcessor.lineEquation(obj.lines);
            
            % Store the laser points
            obj.M = M;
            
            % Create a cell array for storing the delineated tree segments
            obj.trees = cell(1,size(obj.lines,1));
            
        end
        
        function delineateTrees(obj,radTH,distTH)
            %delineateTrees Delineates tree segments around detected lines.
            %   delineateTrees(obj,radTH,distTH) Takes two input arguments:
            %   radTH: The radius threshold. Laser points whose distance to
            %   a line is less than or equal to radTH are automatically
            %   assigned to the line.
            %   distTH: The distance threshold. Laser points whose distance
            %   to an already assigned point is less than distTH are
            %   assigned to the same line as this point.
            %
            %   The function delineates tree segments by assigning laser
            %   points to lines detected earlier from the point cloud. The
            %   function goes through all detected lines in a descending
            %   order based on the length of the lines. All laser points
            %   whose distance to a line is less than radTH are assigned to
            %   this line. Then, using iterative region growing, additional
            %   points are added to the line by comparing the distances
            %   between points already assigned to the line and points that
            %   are not yet assigned to any line. Once a point is assigned
            %   to a line, it is removed from the array of "free" points.
            
            % Go through all lines
            for i = 1:size(obj.lines,1)
                
                % Assign all close points to the line
                pts = obj.findClosePoints(obj.lines(i,:),radTH);
                
                % Assign points to the segment with an iterative region
                % growing process
                pts = obj.regionGrowing(pts,distTH);
                
                obj.trees{i} = pts;
            end

        end
		
		function removeFalseTrees(obj,net,cellSize,cutoffProb)
			% removeFalseTrees Removes the tree segments that are
			% classified as not trees by a classifier
			%	removeFalseTrees(obj,net,cellsSize,cutoffProb) Takes three
			%	input arguments:
			%	net: A trained neural network used for classifying the
			%	segments
			%	cellSize: The cell size of the grids created from the
			%	segments
            %   cutoffProb: The cutoff probability. Segments whose
            %   probability of being a tree is predicted to be larger than
            %   or equal to the cutoff probability are classified as trees.
            %   The smaller the cutoff probability the more segments are
            %   classified as trees. The cutoff probability must be in the
            %   range [0,1].
			%
			%	The function classifies the object's tree segments as tree
			%	or not tree. The points, slopes, lines and intercepts of
			%	the trees classified as not tree are removed from the
			%	object.
            
            % Transform the coordinate system of the tree segments so that
            % the tree segments are always vertical
            transformedTrees = obj.transformSegments;
            
			% Convert the tree segments into raster format
            gridsArray = processSegments(transformedTrees,cellSize,...
                [227,227]);
			
            % Predict whether the segment is a tree or not 
			predictions = predict(net,gridsArray);
            % Use only the probabilites of being a tree, as it contains all
            % the required information (probability of not being a tree is
            % 1-probability of being a tree)
            predictions = predictions(:,2);
            
            % Remove all segments and their corresponding lines, slopes and
            % intercepts whose probability of being a tree is smaller than
            % the cutoff probability
            notTreeLocs = predictions < cutoffProb;
            obj.lines(notTreeLocs,:) = [];
            obj.slopes(notTreeLocs,:) = [];
            obj.intercepts(notTreeLocs,:) = [];
            obj.trees(notTreeLocs) = [];
        end
        
        function transformedSegments = transformSegments(obj)
            %transformSegments For each tree segment, creates a coordinate
            %transformation in which the line representation of the segment
            %lies on the y-axis. Transforms the points belonging to each
            %segment with the corresponding coordinate transformation.
            %   transformedSegments = transformSegments(obj) Takes no input
            %   arguments
            %
            %   The original segments stored in obj.trees are not
            %   transformed. The function returns the transformed points as
            %   a cell array.
            
            transformedSegments = cell(size(obj.trees));
            for t = 1:length(obj.trees)
                % Transform the points of one segment
                transformedSegments{t} = LineProcessor.transformPoints(...
                    obj.lines(t,:),obj.trees{t});
            end
        end
		
        function plotTrees(obj,newFigure,separate,pointSize)
            % plotTrees Plots each delineated tree segment with an
            % individual color.
            %   plotTrees(obj,newFigure,separate,pointSize) Takes three input
            %   arguments:
            %   newFigure: A boolean defining whether the trees should be
            %   plotted to a new figure
            %   separate: A boolean defining whether each tree should be
            %   plotted to its own figure.
            %   pointSize: The plotting size of the points
            
            if separate
                for i = 1:length(obj.trees)
                    figure
                    scatter(obj.trees{i}(:,1),obj.trees{i}(:,2),...
                        pointSize,'filled')
                    axis equal
                end
            else
                if newFigure
                    figure
                end
                
                for i = 1:length(obj.trees)
                    scatter(obj.trees{i}(:,1),obj.trees{i}(:,2),...
                        pointSize,'filled')
                    hold on
                end
                
                axis equal
            end
        end
        
    end
    
    methods (Access = private)
        
        function pts = findClosePoints(obj,line,TH)
            %findClosePoints Finds the points in obj.M that are located
            %close to the given line.
            %   pts = findClosePoints(obj,line,TH) Takes two input arguments:
            %   line: The end points of the line to which the points will
            %   be assigned.
            %   TH: The distance threshold. Points whose distance to the
            %   line is less than TH are determined as close points.
            
            % Find the points that are close to the line
            % segment
            locs = LineProcessor.distPointsSegments(obj.M,line) <= TH;
            pts = obj.M(locs,:);
            
            % Remove the found points from obj.M
            obj.M(locs,:) = [];
        end
        
        function pts = regionGrowing(obj,pts,TH)
            %regionGrowing Assigns "free" points to a segment based on
            %their minimum distance to the segment.
            %   pts = regionGrowing(obj,pts,TH) Takes two input arguments:
            %   pts: An array containing points belonging to a segment
            %   TH: A distance threshold. Points that whose minimum
            %   distance to the segment is smaller than or equal to the
            %   threshold are assigned to the segment.
            %
            %   The function assigns the "free" points to the segment
            %   iteratively. First, it searches for all points whose
            %   minimum distance to the segment is smaller than or equal to
            %   the threshold and assigns these points to the segment.
            %   Then, it repeats this process until no more free points
            %   close to the segment are found.
            
            % Assign points to the segment with an iterative region
            % growing process
            flag = true;
            while flag
                % Calculate the minimum distances between the not
                % assigned and assigned points
                dists = min(pdist2(obj.M(:,1:2),pts(:,1:2)),[],2);

                % Find the indices of the not assigned points that are
                % to be assigned to the segment
                locs = dists <= TH;

                % If some such points were found
                if sum(locs) > 0 
                    % Assign the close points to the segment
                    pts = [pts;obj.M(locs,:)];
                    % Remove the newly assigned points from obj.M
                    obj.M(locs,:) = [];
                % If no close points were found
                else
                    % Move on to the next line segment
                    flag = false;
                end
            end
        end
        
    end
        
end

