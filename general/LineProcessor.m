classdef LineProcessor
    %LineProcessor Contains methods for processing lines and line segments
    
    methods (Static)
       
        % 2D line segment functions
        function slopeDegrees = lineSlope(lineSegments)
            %Calculates the slope of the given line segments. The slope is
            %given as an angle between the line and the x-axis. The slope
            %is in the range [-90,90]

            slopeDegrees = atand( (lineSegments(:,4)-...
                lineSegments(:,3))./...
                (lineSegments(:,2)-lineSegments(:,1)) );
        end
        
        function [slopes,intercepts] = lineEquation(lineSegments)
            %lineEquation Calculates the line equations for the given line
            %segments
            %   [slopes,intercepts] = lineEquation(lineSegments) Takes one
            %   input argument:
            %   lineSegments: An array containing line segments. The format
            %   of the line segments is [startX endX startY endY].
            %
            %   The function returns the slope and intercept of each line
            %   segment. If the line is vertical, the function returns the
            %   slope Inf and the intersection of the line and the x axis.
            
            slopes = ( lineSegments(:,4)-lineSegments(:,3) ) ./...
                ( lineSegments(:,2)-lineSegments(:,1) );
            
            % Ensure that all vertical lines have an Inf slope (and not
            % -Inf)
            slopes(abs(slopes) == Inf) = Inf;
            
            intercepts = lineSegments(:,3) - slopes.*lineSegments(:,1);
            
            % Intercept for vertical lines
            intercepts(abs(intercepts) == Inf) =...
                lineSegments(abs(intercepts) == Inf,1);
        end
        
        function diffs = slopeDiffs(lineSegments,referenceLine)
            %slopeDiffs Calculates the slope differences (in degrees)
            %between the line segments and the reference line
            %   diffs = slopeDiffs(lineSegments,referenceLine) Takes two
            %   input arguments:
            %   lineSegments: An array containing the end points of line
            %   segments
            %   referenceLine: A vector containing the end points of the
            %   reference line.
            
            % Calculate the slopes (in degrees) of the line segments
            segmentSlopes = LineProcessor.lineSlope(lineSegments);
            
            % Calculate the slope (in degrees) of the reference line
            referenceSlope = LineProcessor.lineSlope(referenceLine);
            
            % Calculate the slope difference. The slope difference is
            % either the difference between the slopes or the sum of
            % their absolute complementary angles (90-abs(slope))
            diffs = min(abs(segmentSlopes-referenceSlope),...
                180-abs(segmentSlopes)-abs(referenceSlope));
        end
        
        function dists = distSegmentsSegment(lineSegments,referenceLine)
            %distsSegmentsSegment Calculates the distances of all line
            %segments to the given reference line segment
            %   dists = distSegmentsSegment(lineSegments,referenceLine)
            %   Takes two input arguments:
            %   lineSegments: An array containing end points of one or more
            %   line segments
            %   referenceLine: A vector containing the end points of the
            %   reference line segment
            %
            %   The function calculates the shortest distances between the
            %   given line segments and the reference line segment. The
            %   shortest distance between two line segments AB and CD is
            %   one of the following:
            %   1. The distance between point A and the line segment CD.
            %   2. The distance between point B and the line segment CD.
            %   3. The distance between point C and the line segment AB.
            %   4. The distance between point D and the line segment AB.
            %   5. The intersection of the two line segments, in which case
            %   the distance is 0.
            %
            %   The distances are returned as a vector in which the
            %   distance in position i represents the distance between the
            %   reference line segment and the line segment in position i
            
            % Calculate the potential minimum distances
            dist1 = LineProcessor.distPointsSegments(...
                lineSegments(:,[1 3]),referenceLine);
            dist2 = LineProcessor.distPointsSegments(...
                lineSegments(:,[2 4]),referenceLine);
            dist3 = LineProcessor.distPointsSegments(...
                referenceLine([1 3]),lineSegments);
            dist4 = LineProcessor.distPointsSegments(...
                referenceLine([2 4]),lineSegments);
            
            dists = min([dist1,dist2,dist3,dist4],[],2);
            
            % Determine if the line segments intersect the reference line
            % (http://www.cs.swan.ac.uk/~cssimon/line_intersection.html)
            
            x1 = lineSegments(:,1);
            x2 = lineSegments(:,2);
            y1 = lineSegments(:,3);
            y2 = lineSegments(:,4);
            x3 = referenceLine(1);
            x4 = referenceLine(2);
            y3 = referenceLine(3);
            y4 = referenceLine(4);
            
            ta = ( (y3-y4)*(x1-x3) + (x4-x3)*(y1-y3) ) ./...
                ( (x4-x3)*(y1-y2) - (x1-x2)*(y4-y3) );
            tb = ( (y1-y2).*(x1-x3) + (x2-x1).*(y1-y3) ) ./...
                ( (x4-x3)*(y1-y2) - (x1-x2)*(y4-y3) );
            
            % Assign the distances of those line segments that intersect
            % the reference line segment to zero
            dists(ta >= 0 & ta <= 1 & tb >= 0 & tb <= 1) = 0;
        
        end
        
        function dists = distPointsSegments(pts,segments)
            %distPointsSegments Calculates the distances between the given
            %points and line segments.
            %   dists = distPointsSegments(pts,segments) Takes two input
            %   arguments:
            %   pts: A 2-dimensional vector/nx2 array containing x and y
            %   coordinates of points
            %   segments: A 4-dimensional vector/mx4 array containing the
            %   end points of line segments.
            %   If pts is an array, segments must be a vector and vice
            %   versa. However, both pts and segments can be a vector.
            %
            %   The function calculates the distances between the given
            %   points and line segments.
            %   Theory can be found here: 
            %   http://paulbourke.net/geometry/pointlineplane/
            
            % Check that the input argument requirements are met
            if min(size(pts)) > 1 && min(size(segments)) > 1
                error('Both pts and segments cannot be arrays.')
            end
            
            x1 = segments(:,1);
            y1 = segments(:,3);
            x2 = segments(:,2);
            y2 = segments(:,4);
            x3 = pts(:,1);
            y3 = pts(:,2);
            
            px = x2-x1;
            py = y2-y1;
            norm = px.*px + py.*py;
            
            u = ( (x3-x1).*(x2-x1) + (y3-y1).*(y2-y1) ) ./...
                ( (x2-x1).^2 + (y2-y1).^2 );
            
            % If the intersection of the perpendicular line is outside the
            % line segment
            u(u>1) = 1;
            u(u<0) = 0;
            
            x = x1+u.*px;
            y = y1+u.*py;
            
            dx = x-x3;
            dy = y-y3;
            
            dists = sqrt(dx.*dx + dy.*dy);
            
        end
        
        function len = lineSegmentLength(lineSegment)
            %lineSegmentLength Calculates the length of the given line
            %segment.
            %   len = lineSegmentLength(lineSegment) Takes one input
            %   argument:
            %   lineSegment: The end points of a line segment in format
            %   [startX endX startY endY]
            
            len = LineProcessor.distPts2Pts(lineSegment([1 3]),...
                lineSegment([2 4]));
        end
        
        function lineSegment = createLineSegment(pt1,pt2)
            %createLineSegment Creates a line segment of the given points
            %   lineSegment = createLineSegment(pt1,pt2) Takes two input
            %   arguments:
            %   pt1: One point
            %   pt2: Another point
            %
            %   The function orders the points into a line segment of
            %   format [startX endX startY endY] where startX <= endX
            
            lineSegment = sortrows([pt1;pt2],1);
            lineSegment = [lineSegment(1,1) lineSegment(2,1)...
                lineSegment(1,2) lineSegment(2,2)];
        end
        
        function boolean = isWithin(lineSegment,slope,intercept,...
                circleCenter,radius)
            %isWithin Finds out if line segment is within the given circle.
            %   boolean = isWithin(lineSegment,slope,intercept,circeCenter,
            %   radius)
            %   Takes five input arguments:
            %   lineSegment: A vector containing the end points of a line
            %   segment
            %   slope: The slope of the line segment. Inf for vertical
            %   lines.
            %   intercept: The intercept of the line segment. x-intercept
            %   for vertical lines, y-intercept for others
            %   circleCenter: The center coordinates of the circle
            %   radius: The circle radius
            %
            %   The function finds out whether the given line segment is at
            %   least partially within the given circle. If so, the
            %   function returns 1. Else, the function returns 0.
        
            % Find the intersection coordinates of the line segment (or its
            % extension) and the circle
            [intX,intY] = linecirc(slope,intercept,circleCenter(1),...
                circleCenter(2),radius);
            
            % If no intersections were found
            if isnan(intX(1))
                boolean = 0;
            
            % If intersections were found
            else
                
                % Find out if intersections are located on the range of the
                % line segment. Inspect the x coordinate range, as startX
                % <= endX in the line segment representation
                
                % If at least one intersection is within the line segment
                % x-range
                if sum(intX >= min(lineSegment(1:2)) & intX <=...
                        max(lineSegment(1:2))) > 0
                    
                    % If the line is not vertical, then the intersection is
                    % within the line segment range
                    if slope ~= Inf
                        boolean = 1;
                    % If the line is vertical
                    else
                        % Check the y-range
                        if sum(intY >= min(lineSegment(3:4)) & intY <=...
                                max(lineSegment(3:4))) > 0
                            boolean = 1;
                        else
                            boolean = 0;
                        end
                    end
                    
                % If no intersections are within the line segment range
                else
                    
                    % Check if the line segment is fully within the circle
                    % by calculating the distance between one end point of
                    % the segment and the circle center
                    
                    % If line segment is inside circle
                    if sqrt( (lineSegment(1)-circleCenter(1))^2 +...
                            (lineSegment(3)-circleCenter(2))^2 ) <= radius
                        boolean = 1;
                    
                    % If line segment is outside circle
                    else
                        boolean = 0;
                    end
                end
                  
            end
            
        end
        
        function inRange = isInRange(segments,refSegment)
            %isInRange Calculates which segments are within the range of
            %the reference segment in the principal direction of the
            %reference segment.
            %   inRange = isInRange(segments,refSegment) Takes two input
            %   arguments:
            %   segments: An array of line segments. Each row of the array
            %   represents one line segment. The columns of the array
            %   represent the start x, en x, start y and end y coordinates
            %   of the segments.
            %   refSegment: A vector representing a reference segment to
            %   which the other segments are compared.
            %
            %   The function first determines the principal direction of
            %   the reference segment (the direction to which the rate of
            %   change is faster, x or y). It then finds out which segments
            %   have at least one end point within the range of the
            %   reference segment.
            
            slopeDegrees = LineProcessor.lineSlope(refSegment);
            
            % X direction is principal direction
            if abs(slopeDegrees) <= 45
                i = 1;
                j = 2;
            % Y direction is principal direction
            else
                i = 3;
                j = 4;
            end
                
            minRef = min(refSegment(i:j));
            maxRef = max(refSegment(i:j));

            inRange1 = segments(:,i) >= minRef & segments(:,i) <= maxRef;
            inRange2 = segments(:,j) >= minRef & segments(:,j) <= maxRef;
            inRange = inRange1 | inRange2;
            
        end
        
        function transformedPoints = transformPoints(lineSegment,pts)
            %transformPoints Rotates and translates the given line segment
            %and the given points so that the line segment lies on the
            %y-axis
            %
            %   transformedPoints = transformPoints(lineSegment,pts) Takes
            %   two input arguments: 
            %   lineSegment: A line segment that determines the rotation.
            %   pts: An array in which each row represents one point
            
            % Coordinates of the origin
            x0 = lineSegment(1);
            y0 = lineSegment(3);
            
            % x and y coordinates of the points
            V = pts(:,1:2)';
            
            % Calculate the angle between the line segment and the x-axis 
            alpha = LineProcessor.lineSlope(lineSegment);
            
            % Calculate the rotation angle
            beta = 90-alpha;
            
            % Rotate and translate the points. The translation is done to
            % the counterclockwise direction.
            pts(:,1:2) = ([cosd(beta) -sind(beta);sind(beta)...
                cosd(beta)]*(V-[x0;y0]))';
            transformedPoints = pts;
        end
        
        function lineSegment = sectionWithinCircle(origLineSegment,...
                circleCenter,radius)
            %sectionWithinCircle Determines the part of the tree segment
            %that is within the given circle. Calculates its volume.
            %   lineSegment = sectionWithinCircle(origLineSegment, 
            %   circeCenter,radius)
            %   Takes three input arguments:
            %   origLineSegment: The coordinates of the end points of a
            %   line segment in format [startX endX startY endY]
            %   circleCenter: The coordinates of the center point of a
            %   circle
            %   radius: The radius of the circle
            %  
            %   The function determines the section of the line segment
            %   that is within the given circle and returns it as a new
            %   line segment. If the line segment is fully outside the
            %   circle, the function returns an empty array.
            
            
            % Calculate the slope and intercept of the tree segment's line
            % representation
            [slope,intercept] = LineProcessor.lineEquation(...
                origLineSegment);
            
            % If the line segment is not even partially within the circle
            if ~LineProcessor.isWithin(origLineSegment,slope,intercept,...
                    circleCenter,radius)
                lineSegment = [];

            
            % If the line segment is at least partially within the circle
            else
                % Calculate the distance from the end points of the line
                % segment to the circle center
                distS = LineProcessor.distPts2Pts(...
                    origLineSegment([1 3]),circleCenter);
                distE = LineProcessor.distPts2Pts(...
                    origLineSegment([2 4]),circleCenter);
                
                % If both points are within the circle
                if distS <= radius && distE <= radius
                    % Return the full line segment
                    lineSegment = origLineSegment;
                
                % If both points are outside the circle
                elseif  distS > radius && distE > radius
                    % Calculate the locations of the intersections of the
                    % line segment and the circle. 
                    [intX,intY] = linecirc(slope,intercept,...
                        circleCenter(1),circleCenter(2),radius);
                    
                    % The intersectiond are the end points of the new line
                    % segment
                    lineSegment = [intX intY];    
                
                % If the start point is inside the circle
                elseif distS <= radius
                    % Determine the intersection of the line segment and
                    % the circle
                    intPt = LineProcessor.calculateIntersection(slope,...
                        intercept,circleCenter(1),circleCenter(2),...
                        radius,origLineSegment([2 4]));
              
                    % The start point of the original line segment and the
                    % intersection are the end points of the new line
                    % segment
                    lineSegment = [origLineSegment(1) intPt(1)...
                        origLineSegment(3) intPt(2)];
                        
                % If the end point is inside the circle
                else
                    % Determine the intersection of the line segment and
                    % the circle
                    intPt = LineProcessor.calculateIntersection(slope,...
                        intercept,circleCenter(1),circleCenter(2),...
                        radius,origLineSegment([1 3]));
              
                    % The end point of the original line segment and the
                    % intersection are the end points of the new line
                    % segment
                    lineSegment = [origLineSegment(2) intPt(1)...
                        origLineSegment(4) intPt(2)];
                    
                end
                
                
                % Order the end points so that the end point with the
                    % smaller x coordinate comes first
                if lineSegment(2) < lineSegment(1)
                    lineSegment = [lineSegment(2) lineSegment(1)...
                        lineSegment(4) lineSegment(3)];
                end
            end
         
        end
        
        
        % Line functions
        function intPt = calculateIntersection(slope,intercept,...
                centerX,centerY,radius,outsideLoc)
            %calculateIntersection Calculates the intersection
            %coordinates of the given line segment and the given sample
            %plot.
            %   intPt = calculateIntersection(slope,intercept,centerX,
            %   centerY,radius,outsideLoc)
            %   Takes six input arguments:
            %   slope: The slope of the line segment
            %   intercept: The intercept of the line segment
            %   centerX: The x coordinate of the center point of the circle
            %   centerY: The y coordinate of the center point of the circle
            %   radius = The radius of the circle
            %   outsideLoc: The end point of the line segment that is
            %   outside the circle
            
            % Calculate the intersection points
            [intX,intY] = linecirc(slope,intercept,centerX,centerY,radius);

            % Find out which intersection point is closer to the outside
            % end of the tree. This location is the correct intersection
            % between the sample plot border and line segment
            dist1 = LineProcessor.distPts2Pts(outsideLoc,...
                [intX(1),intY(1)]);
            dist2 = LineProcessor.distPts2Pts(outsideLoc,...
                [intX(2),intY(2)]);

            if dist1 < dist2
                intPt = [intX(1),intY(1)];
            else
                intPt = [intX(2),intY(2)];
            end

        end 
        
        function projectedPts = projectToLine(pts,A,B)
            %Projects points to the given line. Returns the coordinates of
            %the projected points.
            %   projectedPts = projectToLine(pts,A,B) Takes three input
            %   arguments:
            %   pts: An array of points. The first two columns of the array
            %   contain the x and y coordinates of the points.
            %   A: The line slope
            %   B: The line intercept
            %
            %   The function projects the input points to the input line
            %   and returns the coordinates of the projected points. Does
            %   not return the information on the additional columns of the
            %   pts array. However, the corresponding points in pts and the
            %   output points are on the same row, so the additional
            %   information can be easily added to the points later on if
            %   needed.
            
            % If the line is not vertical or horizontal
            if abs(A) ~= Inf && A ~= 0
                projectedPts(:,1) = (pts(:,2)+1/A*pts(:,1)-B)/(A+1/A);
                projectedPts(:,2) = (B+A*(pts(:,2)+1/A*pts(:,1)-B)/...
                    (A+1/A));
            % If the line is vertical
            elseif abs(A) == Inf
                projectedPts(:,1) = ones(size(pts,1),1)*B;
                projectedPts(:,2) = pts(:,2);
            % If the line is horizontal    
            else
                projectedPts(:,1) = pts(:,1);
                projectedPts(:,2) = ones(size(pts,1),1)*B;
            end
            
        end
        
        function slopeDegrees = slopeToAngle2D(slope)
            %slopeToAngle2D Calculates the rotation angle of the given
            %slope
            %   slopeDegrees = slopeToAngle2D(slope) Takes one input
            %   argument:
            %   slope: The slope of a line
            %   
            %   The function converts the given slope to the corresponding
            %   rotation angle. The angle is given as the angle between the
            %   positive x axis and the given line. The angle is between
            %   [-90,90].
            
            slopeDegrees = atand(slope);
        end
        
        
        % Point functions
        function dist = distPts2Pts(pt1,pt2)
            %distPts2Pts Calculates the row-wise distances between points
            %   dist = distPts2Pts(pt1,pt2) Takes two input arguments:
            %   pt1: An  nx2 array in which each row contains the x and y
            %   coordinate of a point
            %   pt2: An nx2 array in which each row contains the x and y
            %   coordinate of a point
            %
            %   The function calculates the row-wise distances between the
            %   points in pt1 and pt2. For example, the distance in
            %   position 3 of the output array represents the distance
            %   between pt1(3,:) and pt2(3,:).
            
            dist = sqrt((pt1(:,1)-pt2(:,1)).^2 + (pt1(:,2)-pt2(:,2)).^2);
        end
        
        % Plotting
        function plotLines(lines,newFig,lineWidth,color)
            %plotLines Plots the given lines
            %   fig = plotLines(lines,newFig,lineWidth,color) Takes
            %   four input arguments:
            %   lines: An array containing the start and end points of
            %   lines to be plotted. Each row of the array contains one
            %   line. The line representation is of the form [startX, endX,
            %   startY, endY].
            %   newFig: A logical value determining whether the lines will
            %   be plotted on a new figure
            %   lineWidth: The width of the lines
            %   color: A string defining the color of the lines. OPTIONAL
            %
            %   The function plots the given lines and returns the plot
            %   handle.
            
            if newFig
                figure;
                hold on
                axis equal
            end
            
            % If the color argument was not given
            if nargin == 3
                for i = 1:size(lines,1)
                    % Plot each line with and individual color
                    plot(lines(i,1:2),lines(i,3:4),'LineWidth',lineWidth)
                end
            else
                for i = 1:size(lines,1)
                    % Plot each line with the same color
                    plot(lines(i,1:2),lines(i,3:4),color,...
                        'LineWidth',lineWidth)
                end
            end
        end
        
        % Writing to file
        function writeShape(lines,filename)
            %writeShape Writes the given lines into a shapefile.
            %   writeShape(lines,filename) Takes two input arguments:
            %   lines: An array containing line segments. Each row of the
            %   array represents one line segment. The line format is
            %   [startX endX startY endY].
            %   filename: The name of the shapefile to be written
            
            if isempty(lines)
                warning(strcat("The input argument lines was empty.",...
                    " No shapefile was written."))
                
            else
                for i = 1:size(lines,1)

                    s(i).ID = i;
                    s(i).Geometry = 'line';
                    s(i).X = lines(i,1:2);
                    s(i).Y = lines(i,3:4);

                end

                S = mapshape(s);
                spec = makedbfspec(S);

                shapewrite(S,filename,'DbfSpec',spec)
            end
        end
        
    end
    
        
end

