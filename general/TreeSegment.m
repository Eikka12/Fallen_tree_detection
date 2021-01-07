classdef TreeSegment
    %TreeSegment A class storing information of a single tree segment
    
    properties
        lineSegment;
        origPts; % Points with their original coordinates
        pts; % Points transformed to lie on the y-axis
        len; % The length of the segment
        diameter; % The diameter of the segment
        
    end
    
    methods
        function obj = TreeSegment(pts,lineSegment)
            %TreeSegment Construct an instance of this class
            %   obj = TreeSegment(pts,lineSegment) Takes two input
            %   arguments:
            %   pts: An array of points belonging to the segment.
            %   lineSegment: A 4-element vector containing the endpoints of
            %   the line representation of the segment
            
            % Store the line segment
            obj.lineSegment = lineSegment;
            % Store the points
            obj.origPts = pts;
            % Create a transformation that rotates the line segment so that
            % it lies on the y-axis. Use the same transformation for the
            % points and store the transformed points.
            obj.pts = LineProcessor.transformPoints(lineSegment,pts);
            
            % Calculate the segment length and diameter
            obj.len = obj.calculateLength;
            obj.diameter = obj.calculateDiameter;
        end
        
        function [lineSegment,volume] = sectionWithinCircle(obj,...
                circleCenter,radius)
            %sectionWithinCircle Determines the part of the tree segment
            %that is within the given circle. Calculates its volume.
            %   [lineSegment,volume] = sectionWithinCircle(obj,circeCenter,
            %   radius)
            %   Takes two input arguments:
            %   circleCenter: The coordinates of the center point of a
            %   circle
            %   radius: The radius of the circle
            %  
            %   The function determines the section of the tree segment
            %   that is within the given circle. The function returns the
            %   section as a line segment and an estimated volume. The
            %   volume as follows: The section of the tree segment that is
            %   located within the circle is divided into parts with height
            %   0.1 m. Then, each part's upper and lower diameter
            %   (Dup,Dlow) is calculated and the volume of the part is
            %   calculated as the volume of a truncated cone with diameters
            %   Dup and Dlow and a height 0.1 m. Finally, the volumes of
            %   each part are added to get the total volume.
            
            % Get the section of the line segment that is within the circle
            lineSegment = LineProcessor.sectionWithinCircle(...
                obj.lineSegment,circleCenter,radius);
            
            % If the line segment is fully outside the circle
            if isempty(lineSegment)
                % The volume is zero
                volume = 0;
            % If the line segment is at least partially within the circle    
            else
                
                % The line segment section is now in the original
                % coordinate system. The section must be transformed into
                % the coordinate system in which obj.lineSegment defines
                % the y-axis, as the points are in this coordinate system
                endpoints = LineProcessor.transformPoints(...
                    obj.lineSegment,[lineSegment([1,3]);...
                    lineSegment([2 4])]);
                % Calculate the diameters of the segment within the circle
                % with 0.1 meter intervals
                volume = obj.calculateVolumeOfSection(...
                    sort(endpoints(:,2)),0.1,1);
            end
        end
        
        function plotSegment(obj)
            % Plots the segment as a scatterplot
            
            figure
            scatter(obj.pts(:,1),obj.pts(:,2),10,'filled')
            axis equal;
        end
        
        function vol = calculateVolume(obj,stepSize,maxDiam)
            %calculateVolume Calculates the volume of the tree segment
            %   vol = calculateVolume(obj,stepSize,maxDiam) Takes two input
            %   arguments:
            %   stepSize: The volume is calculated by dividing the tree
            %   segment into intervals with height stepSize (in meters).
            %   maxDiam: The maximum diameter. If the diameter at interval
            %   HI is over the maximum diameter, this diameter is ignored
            %   and the parts of the tree segment this height interval
            %   divides are merged in the volume calculations. THE
            %   PARAMETER IS OPTIONAL. The default value for maxDiam is 1.
            %
            %   The function calculates the volume of the tree segment by
            %   computing the boundary of the segment, calculating the
            %   boundary diameter at several locations and calculating the
            %   volume as piecewise cut cone volumes
            
            if isempty(obj.pts)
                vol = 0;
                return;
            end
            
            if nargin == 2
                maxDiam = 1;
            end 
            
            % Calculate the volume
            vol = obj.calculateVolumeOfSection([min(obj.pts(:,2))...
                max(obj.pts(:,2))],stepSize,maxDiam);
        end
        
        function vol = calculateVolumeOfSection(obj,section,stepSize,...
                maxDiam)
            %calculateVolumeOfSection Calculates the volume of the defined
            %section of the tree segment using a piecewise truncated cone
            %approach
            %   vol = calculateVolumeOfSection(obj,section,stepSize, 
            %   maxDiam)
            %   Takes three input arguments:
            %   section: A 2-element vector containing the minimum and
            %   maximum y coordinates of the section
            %   stepSize: The height of each truncated cone. The section is
            %   divided into range(section)/stepSize parts for which the
            %   volume is calculated separately. However, the actual height
            %   of each truncated cone may vary depending of the maxDiam
            %   parameter
            %   maxDiam: The maximum allowed diameter of the tree segment.
            %   If the diameter calculated at some point of the tree
            %   segment is over maxDiam, this diameter is ignored and the
            %   two parts of the tree segment which this diameter divides
            %   are merged.
            
            if isempty(obj.pts)
                vol = 0;
                return;
            end
            
            % Calculate the diameters of the segment within the circle
            % with 0.1 meter intervals
            diamHeights = [section(1):stepSize:section(2),section(2)];

            % Calculate the diameters at each height interval and compute
            % the volume between consecutive height intervals (except the
            % last two height intervals). Add them to the total volume
            diams = zeros(size(diamHeights));
            for HI = 1:length(diamHeights)
                % Calculate the diameter at the current height interval
                diams(HI) = obj.calculateDiameterAtHeight(...
                    diamHeights(HI));
            end

            % Find all diameters that are over the maximum diameter
            locs = diams > maxDiam;

            % Remove these diameters and corresponding height intervals
            % from the volume calculations
            diams(locs) = [];
            diamHeights(locs) = [];

            % Calculate the total volume as the sum of volumes of the
            % sections between each height interval
            vol = 0;
            for HI = 2:length(diamHeights)
                vol = vol + obj.cutConeVol(diams(HI-1),diams(HI),...
                    diamHeights(HI)-diamHeights(HI-1)); 
            end
        end
        
        function diam = calculateDiameterAtMidpoint(obj)
            % calculateDiameterAtMidpoint Calculates the tree segment
            % diameter at the midpoint of the tree
            %   diam = calculateDiameterAtMidpoint(obj) Takes no input
            %   arguments
            %
            %   The function calculates the midpoint diameter as the
            %   average diameter of a one-meter section of the midpoint of
            %   the tree sampled at 0.01 meter intervals
            
            if isempty(obj.pts)
                diam = 0;
                return;
            end
            
            % Calculate the midpoint of the tree
            midpoint = (max(obj.pts(:,2))-min(obj.pts(:,2)))/2;
            
            % Calculate the diameters of the one-meter section around the
            % midpoint
            diamPoints = (midpoint-0.5):0.01:(midpoint+0.5);
            diams = zeros(size(diamPoints));
            for i = 1:length(diamPoints)
                diams(i) = obj.calculateDiameterAtHeight(diamPoints(i));
            end
            
            % Calculate the midpoint diameter as the average of the
            % previously calculated diameters
            diam = mean(diams);
        end
                
        function diam = calculateDiameterAtHeight(obj,y)
            % calculateDiameterAtHeight Calculates the diameter at the
            % specified y coordinate location
            %   diam = calculateDiameterAtHeight(obj,y) Takes one
            %   input argument:
            %   y: The y coordinate at which the diameter is calculated.
            %
            %   The function calculates the diameter at the specified
            %   height by computing the boundary of the segment, drawing a
            %   horizontal line at the specified height and computing the
            %   distance between the left and right intersection points of
            %   the horizontal line and the boundary.
            
            % Compute the boundary of the tree segment and extract the
            % boundary points
            bdaryPts = obj.pts(boundary(obj.pts(:,1:2)),1:2);
            
            % Get the intersection facets of the boundary and the
            % horizontal line
            [facet1,facet2] = obj.findIntersectionFacets(bdaryPts,y);
            
            % If there are not intersection facets
            if isempty(facet1) || isempty(facet2)
                diam = 0;
            else
                % Calculate the x-coordinates of the intersection points of
                % the facets and the horizontal line
                [slopes,intercepts] = LineProcessor.lineEquation(...
                    [facet1;facet2]);
                
                intersectionXs = (y-intercepts)./slopes;
                
                % If either of the line segments is vertical, the
                % intersection's x coordinate is the intercept of this line
                % segment (LineProcessor returns the x intercept in the
                % case of vertical lines
                nanLocs = abs(slopes) == Inf;
                intersectionXs(nanLocs) = intercepts(nanLocs);
                
                % The diameter is the distance between the intersection
                % points
                diam = abs(intersectionXs(1)-intersectionXs(2));
                
                % However, if either of the intersection facets is
                % horizontal, the diameter is 0
            end
        end
        
        function len = calculateLength(obj)
            %calculateLength Calculates the length of the segment
            %   len = calculateLength(obj) Takes no input arguments
            %
            %   The function calculates the length of the tree segment as
            %   the range of the points beloging to the segment in the
            %   y-axis direction. NOTE THAT THE SEGMENT HAS BEEN COORDINATE
            %   TRANSFORMED IN A WAY THAT THE TREE SEGMENT IS ALIGNED WITH
            %   THE Y-AXIS
                
            if size(obj.pts,1) > 1    
                len = range(obj.pts(:,2));
            else
                len = 0;
            end
        end
                    
    end
    
    methods (Access = private)
        
        function diam = calculateDiameter(obj)
            %calculateDiameter Calculates the diameter of the segment.
            %   diam = calculateDiameter(obj) Takes no input arguments.
            %
            %   The function calculates the diameter of the tree segment.
            %   First, the function extracts the points that are located
            %   more than 1.5 meters away (in the y direction) from the
            %   minimum and maximum y coordinates of the tree. The purpose
            %   of this step is to remove the possibly included tree stump
            %   that might distort the further measurements. The point
            %   extraction step is performed only if the tree length is
            %   more than 4 meters. Next, the function determines the area
            %   of the boundary polygon of the remaining laser points.
            %   Finally, the tree diameter is determined as the width of a
            %   rectangle with an area equal to the area of the boundary
            %   polygon and height equal to the range of the laser points
            %   in y direction. NOTE THAT THE SEGMENT HAS BEEN COORDINATE
            %   TRANSFORMED IN A WAY THAT THE TREE SEGMENT IS ALIGNED WITH
            %   THE Y-AXIS
            
            points = obj.pts;
            
            % If the segment contains more than two points
            if size(points,1) > 2
                % Get the maximum and minimum y coordinates of the tree
                % points
                maxY = max(points(:,2));
                minY = min(points(:,2));

                % If the length of the tree is more than 4 meters
                if obj.calculateLength > 4

                    % Get the tree points that are more than 1.5 m away
                    % from the tree top and tree stump in the y direction
                    locs = points(:,2) < maxY - 1.5 & points(:,2) >...
                        minY + 1.5;
                    points = points(locs,:);
                end

                % Determine the boundary area of the boundary polygon of
                % the remaining points
                [~,A] = boundary(points(:,1:2));

                % Calculate the diameter as the width of a rectangle with
                % the same area as the area of the boundary polygon and a
                % length equal to the range of the remaining points in y
                % direction
                diam = A/range(points(:,2));
            else
                diam = 0;
            end
        end
        
    end
    
    methods (Static,Access = private)
        
        function [facet1,facet2] = findIntersectionFacets(bdaryPts,y)
            %findIntersectionFacets Finds the boundary facets that
            %intersect the horizontal line at y.
            %   [facet1,facet2] = findIntersectionFacets(bdaryPts,y) Takes
            %   two input arguments:
            %   bdaryPts: An array of points representing the boundary of
            %   the tree segment
            %   y: The y coordinate of the horizontal line.
            %
            %   The function searches for the two facets of the boundary of
            %   the tree segment that intersect with the horizontal line at
            %   y.
            
            facet1 = [];
            facet2 = [];
            for pt = 2:size(bdaryPts,1)
                % If the two consecutive boundary points are on the two
                % sides of the horizontal line
                if bdaryPts(pt-1,2) >= y && bdaryPts(pt,2) <= y
                    % If the found boundary points have the same y
                    % coordinate.
                    if bdaryPts(pt-1,2) == bdaryPts(pt,2)
                        % Continue to the next iteration without storing
                        % the facet
                        continue;
                    else
                        % If this is the first intersection facet found
                        if isempty(facet1)
                            % Store the facet
                            facet1 = LineProcessor.createLineSegment(...
                                bdaryPts(pt-1,1:2),bdaryPts(pt,1:2)); 
                        % If this is the second intersection facet found
                        else
                            % Store the facet
                            facet2 = LineProcessor.createLineSegment(...
                                bdaryPts(pt-1,1:2),bdaryPts(pt,1:2)); 
                            % End processing
                            break;
                        end
                    end
                elseif bdaryPts(pt-1,2) <= y && bdaryPts(pt,2) >= y
                    % If this is the first intersection facet found
                    if isempty(facet1)
                        % Store the facet
                        facet1 = LineProcessor.createLineSegment(...
                            bdaryPts(pt-1,1:2),bdaryPts(pt,1:2)); 
                    % If this is the second intersection facet found
                    else
                        % Store the facet
                        facet2 = LineProcessor.createLineSegment(...
                            bdaryPts(pt-1,1:2),bdaryPts(pt,1:2)); 
                        % End processing
                        break;
                    end
                end
            end
        end
        
        function vol = cutConeVol(diam1,diam2,h)
            %cutConeVol Calculates the volume of a cut cone
            %   vol = cutConeVol(diam1,diam2,h) Takes three input
            %   arguments:
            %   diam1,diam2: The diameters of the cut cone
            %   h: The height of the cut cone
            
            vol = 1/3*pi*h*( (diam1/2)^2 + diam1*diam2/4 + (diam2/2)^2 );
        end
        
    end
end

