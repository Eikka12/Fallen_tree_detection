classdef LineSearcher
    %LineSearcher Contains methods for searching lines and line segments
    %from a point cloud as well as methods for plotting and storing the
    %lines into files.
    
    methods
        
        function obj = LineSearcher()
            %LineSearcher Construct an instance of this class
        end
        
    end
    
    methods (Static)
        
        % New version of the houghLines function. After an infinite line is
        % found, the function forms a line segment to the location that
        % contains the most points close to the line.
        function lines = houghLines(pts,angleInterval,numOffsets,...
                minPeakVal,pDist,maxSeparation)
            %houghLines Finds lines from the given point cloud using the
            %Hough transformation.
            %   lines = houghLines(pts,angleInterval,numOffsets,
            %   minPeakVal,pDist,maxSeparation)
            %   Takes six input arguments:
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
            %   The function finds lines from the given point cloud using
            %   an iterative Hough transformation. The Hough transformation
            %   determines line parameters as peaks in the binned parameter
            %   space. See function houghpoint for more information. At
            %   each iteration, the function searches for the strongest
            %   peak of the Hough transformation and extracts the
            %   parameters of the corresponding line. Then it searches for
            %   the points that are close to this line and removes them
            %   from the point cloud. Once no strong peaks are found, the
            %   function stops the line search. After this, the found
            %   infinite lines are split into line segments based on the
            %   locations of the points that are close to each line.
            %
            %   The line segments are returned as an array in which each
            %   row represents one line. The line representation is of the
            %   form [startX endX startY endY].

            % Create an array for storing the line parameters
            lines = [];

            % Convert the angle interval from degrees to radians
            angleInterval = deg2rad(angleInterval);

            flag = true;
            while flag
                
                % Normalize the data
                ptsNorm = pts;
                minPtsX = min(pts(:,1));
                minPtsY = min(pts(:,2));
                ptsNorm(:,1) = pts(:,1) - minPtsX;
                ptsNorm(:,2) = pts(:,2) - minPtsY;

                % Calculate Hough transformation
                [H,angles,offsets] = houghpoint(ptsNorm(:,1:2),...
                    -pi/2:angleInterval:pi/2,numOffsets);

                % If line search should continue
                if max(H(:)) >= minPeakVal
                    % Get the location of the largest peak of the Hough
                    % transformation
                    [i,j] = find(H == max(H(:)));
                    i = i(1);
                    j = j(1);

                    % Get the slope and intercept of the peak
                    angle = angles(i);
                    offset = offsets(j);
                    
                    % If the line is not vertical
                    if abs(angle) < 1.570
                        % slope
                        A = tan(angle);
                        % intersect + denormalization
                        B = offset/cos(angle) + minPtsY - A*minPtsX;
                    % If the line is vertical
                    else
                        % slope
                        A = angle/1.5708 * Inf;
                        % location where line cuts x axis
                        B = minPtsX + abs(offset);
                    end

                    % Find the points that are close to the line
                    dists = LineSearcher.distPtsLine(pts,A,B);
                    linePts = pts(dists <= pDist,:);
                    
                    % Project the close points to the line
                    projectedPts = LineSearcher.projectToLine(linePts,A,B);
                    
                    % Sort the points in ascending order based on the x
                    % coordinate (or y coordinate if the line is vertical)
                    
                    projectedPts = sortrows(projectedPts);

                    % Calculate the separation between each point going
                    % from left to right
                    consecutiveDists = sqrt(sum((projectedPts(2:end,:)-...
                    projectedPts(1:end-1,:)).^2,2));
                
                    % Find the cut locations
                    cutInds = find(consecutiveDists > maxSeparation);
                    
                    % If there are no cuts in the line
                    if isempty(cutInds)
                        % Create a line segment that spans through all
                        % points
                        segment = [projectedPts(1,1),projectedPts(end,1),...
                            projectedPts(1,2),projectedPts(end,2)];
                        
                        % Store the line segment and delete the points that
                        % are close to it
                        [pts,lines] = LineSearcher.storeAndDelete(pts,A,...
                            B,segment,pDist,lines);
                        
                    % If there is one cut in the line
                    elseif length(cutInds) == 1
                        % Find the largest points group of the line points
                        % and create a line segment from it
                        
                        % If the first group is larger
                        if cutInds > size(projectedPts,1) - cutInds
                            segment = [projectedPts(1,1),...
                                projectedPts(cutInds,1),...
                                projectedPts(1,2),projectedPts(cutInds,2)];
                            
                            % Store the line segment and delete the points
                            % that are close to it
                            [pts,lines] =...
                                LineSearcher.storeAndDelete(pts,A,B,...
                                segment,pDist,lines);
                            
                        % If the second group is larger
                        elseif cutInds < size(projectedPts,1) - cutInds
                            segment = [projectedPts(cutInds+1,1),...
                                projectedPts(end,1),...
                                projectedPts(cutInds+1,2),...
                                projectedPts(end,2)];
                            
                            % Store the line segment and delete the points
                            % that are close to it
                            [pts,lines] =...
                                LineSearcher.storeAndDelete(pts,A,B,...
                                segment,pDist,lines);
                            
                        % If the groups are equally large
                        else
                            % Create a line segment from both groups
                            
                            % First segment
                            segment = [projectedPts(1,1),...
                                projectedPts(cutInds,1),...
                                projectedPts(1,2),projectedPts(cutInds,2)];
                            
                            % Store the line segment and delete the points
                            % that are close to it
                            [pts,lines] =...
                                LineSearcher.storeAndDelete(pts,A,B,...
                                segment,pDist,lines);
                          
                            % Second segment
                            segment = [projectedPts(cutInds+1,1),...
                                projectedPts(end,1),...
                                projectedPts(cutInds+1,2),...
                                projectedPts(end,2)];
                            
                            % Store the line segment and delete the points
                            % that are close to it
                            [pts,lines] =...
                                LineSearcher.storeAndDelete(pts,A,B,...
                                segment,pDist,lines);
                        end
                        
                        
                    % If there are many cuts in the line
                    else
                        % Find the running number of the largest group of
                        % the line points
                        groupSizes = [cutInds;size(projectedPts,1)]-...
                            [0;cutInds];
                        largestRunningNumbers =...
                            find(groupSizes == max(groupSizes));
                        
                        % Include the last points index (number of points)
                        % in the cutInds vector
                        cutInds = [cutInds;size(projectedPts,1)];
                        
                        % Go through all groups with the largest number of
                        % points
                        for i = 1:length(largestRunningNumbers)
                            
                            % If first group is largest group
                            if largestRunningNumbers(i) == 1
                                s = 1;
                                e = cutInds(largestRunningNumbers(i));
                            else
                                s = cutInds(largestRunningNumbers(i)-1)+1;
                                e = cutInds(largestRunningNumbers(i));
                            end
                            
                            % Create a line segment
                            segment = [projectedPts(s,1),...
                                projectedPts(e,1),projectedPts(s,2),...
                                projectedPts(e,2)];
                            
                            % Store the line segment and delete the points
                            % that are close to it
                            [pts,lines] =...
                                LineSearcher.storeAndDelete(pts,A,B,...
                                segment,pDist,lines);
                            
                        end
                            
                    end
                    
                    % If there are no points left to process
                    if isempty(pts)
                        flag = false;
                    end
                    
                % If line search should end   
                else
                    flag = false;
                end
            end
        end
        
        function cleanedLines = lineCleaning(lines,shortTH,slopeDiff,...
                distDiff,overlapTH)
            %lineCleaning Merges and deletes the given lines based on
            %certain criteria.
            %   cleanedLines = lineCleaning(lines,shortTH,slopeDiff,
            %   distDiff,overlapTH)
            %   Takes five input arguments:
            %   lines: A cell array in which each cell stores the lines of
            %   one cell of the gridded point cloud. The lines are stored
            %   as an array in which each row represents one line. The line
            %   format is [startX endX startY endY].
            %   shortTH: A threshold line length. If the length of a line
            %   is shorter than or equal to the threshold, the line is a
            %   short line. Otherwise, the line is a long line.
            %   slopeDiff: The maximum slope difference for the line
            %   segments to be considered belonging to the same line.
            %   distDiff: The maximum distance from a start or end point of
            %   one line segment to the start or end point of another line
            %   segment for the lines to be considered belonging to the
            %   same line.
            %   overlapTH: The overlap threshold (in percentages) used for
            %   determining whether two line segments can be merged. Lines
            %   with x or y range overlaps over the threshold will not be
            %   merged.
            %
            %   The function merges and removes lines from the given lines
            %   array. Line segments whose:
            %   1. slopes differ less than slopeDiff from each other
            %   2. distance differs less than distDiff
            %   3. range overlap percentage is less than overlapTH
            %   will be merged. After merging, the remaining lines that are
            %   shorter than or equal to shortTH will be removed.
            %
            %   The function returned the resulting lines as an array. The
            %   line format is the same as in the input array.

            % Create an array for storing the new lines
            cleanedLines = [];

            % Get all lines to a single array and sort them in ascending
            % order based on the start X coordinate
            lines = LineSearcher.combineLineList(lines);
            if isempty(lines)
                return
            end
            
            lines = sortrows(lines,1);

            % Calculate the line lengths and add them as the last column
            lines(:,5) = LineSearcher.lineLength(lines);
            % Remove lines that have a zero length
            lines(lines(:,5) == 0,:) = [];
            % Calculate the line slopes and add them as the last column
            lines(:,6) = LineSearcher.lineSlope(lines);

            while size(lines,1) >= 1  
                % Find the lines that can potentially be merged with
                % the current line
                [locsPotential,slopeDiffs,dists] =...
                    LineSearcher.findPotentialMergeables(lines,...
                    slopeDiff,distDiff,overlapTH);

                % If no potentially mergeable lines were found
                if sum(locsPotential) == 0 
                    % If the line is a long line, include it in the
                    % output lines
                    if lines(1,5) > shortTH
                        cleanedLines = [cleanedLines;lines(1,:)];
                    end

                    % Remove the line from the lines array
                    lines(1,:) = [];
                
                % If one potentially mergeable line was found
                elseif sum(locsPotential) == 1
                    % Merge the line with the current line, set the merged
                    % line as the first line of the lines array and remove
                    % the original lines from the array.
                    lines = LineSearcher.mergeAndRemove(lines,...
                        locsPotential);
                    
                % If more than one potentially mergeable lines were
                % found
                else
                    
                    % Find the best match of the potentially mergeable
                    % lines
                    cost = (slopeDiffs+1).*dists;
                    % Assure that a not potential line isn't selected as
                    % the best match line
                    cost(~locsPotential) = Inf;
                    bestLoc = cost == min(cost);
                    
                    % If there are several best match lines, select the
                    % first one
                    if sum(bestLoc) > 1
                        Is = find(bestLoc == 1);
                        bestLoc(Is(2:end)) = 0;
                    end
                    
                    % Merge the line with the current line, set the merged
                    % line as the first line of the lines array and remove
                    % the original lines from the array.
                    lines = LineSearcher.mergeAndRemove(lines,bestLoc);
                    
                end

            end
            
            if ~isempty(cleanedLines)
                % Remove short lines from the new lines
                cleanedLines(cleanedLines(:,5) <= shortTH,:) = [];
                
                % Remove the line length and slope from the new lines array
                cleanedLines = cleanedLines(:,1:4);
            end
        end
        
        function combinedLines = combineLineList(lines)
            %Takes the arrays from each cell of the cell array and joins
            %them to a single array in which each row represents the end
            %points of one line.

            combinedLines = [];

            for i=1:length(lines)
                combinedLines = [combinedLines;lines{i}];
            end

        end
        
    end
    
    methods (Static, Access = private)
        
        function [pts,deletedPts] = deleteLineSegmentPoints(pts,A,B,...
                segment,maxDist)
            %deleteLineSegmentPoints Removes the points that are close to
            %the input line segment.
            %   [pts,deletedPts] = deleteLineSegmentPoints(pts,A,B,
            %   segment,maxDist)
            %   Takes five input arguments:
            %   pts: An array of points. The first to columns of the array
            %   contain the x and y coordinates of the points.
            %   A: The slope of the line segment.
            %   B: The intersect of the line segment.
            %   segment: The end points of the line segment. The format of
            %   the segment is [startX endX startY endY].
            %   maxDist: The maximum distance of a point to the line for
            %   the point to be considered as a close point.
            %
            %   The function searches for the points that are close to the
            %   given line segment.
            %
            %   The function returns a points array from which the close
            %   points have been removed.
            
            % Get the points' distances to the line
            dists = LineSearcher.distPtsLine(pts,A,B);
            
            % Find the points that are close to the line
            closeLineInds = dists <= maxDist;
            
            % Find the points that are within the segment range
            
            % If the line is not vertical
            if abs(A) ~= Inf
                
                projectedPts = LineSearcher.projectToLine(pts,A,B);
                withinRangeInds = projectedPts(:,1) >= segment(1) &...
                    projectedPts(:,1) <= segment(2);
            % If the line is vertical
            else
                withinRangeInds = pts(:,2) >= segment(3) & pts(:,2)...
                    <= segment(4);
            end
            
            % Find the indices of the points that are close to the line
            % segment
            closeInds = closeLineInds & withinRangeInds;
            
            % Delete the points that are close to the line segment
            deletedPts = pts(closeInds,:);
            pts = pts(~closeInds,:);
        end
        
        function [pts,lines] = storeAndDelete(pts,A,B,segment,maxDist,...
                lines)
            % Stores the given line segment and deletes points that are
            % within the defined distance from it.
            
            % Delete the points that are close to the line
            [pts,~] = LineSearcher.deleteLineSegmentPoints(pts,A,B,...
                segment,maxDist);

            % Add the segment to the lines array
            lines = [lines; segment];
        end
        
        function dists = distPtsLine(pts,A,B)
            %Calculates the distances from the given points to the given
            %line.
            %   dists = distPtsLine(pts,A,B)
            
            % If the line is not vertical or horizontal
            if abs(A) ~= Inf && A ~= 0
                dists = sqrt((pts(:,1)-((pts(:,2)+1/A*pts(:,1)-B)/...
                    (A+1/A))).^2 + (pts(:,2)-(B+A*(pts(:,2)+1/...
                    A*pts(:,1)-B)/(A+1/A))).^2);
            % If the line is vertical
            elseif abs(A) == Inf
                dists = abs(B-pts(:,1));
            % If the line is horizontal
            else
                dists = abs(B-pts(:,2));
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
            %   B: The line intersect
            %
            %   The function projects the input points to the input line
            %   and returns the coordinates of the projected points.
            
            % If the line is not vertical or horizontal
            if abs(A) ~= Inf && A ~= 0
                projectedPts(:,1) = (pts(:,2)+1/A*pts(:,1)-B)/(A+1/A);
                projectedPts(:,2) = (B+A*(pts(:,2)+1/A*pts(:,1)-B)/(A+1/A));
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
        
        function lengths = lineLength(lines)
            %Calculates the length of the given lines

            lengths = sqrt((lines(:,1)-lines(:,2)).^2 +...
                (lines(:,3)-lines(:,4)).^2);
        end

        function slopes = lineSlope(lines)
            %Calculates the slope of the given lines. The slope is given as
            %an angle between the line and the x-axis.

            slopes = atand( (lines(:,4) - lines(:,3))./...
                (lines(:,2)-lines(:,1)) );
        end
        
        function [locsPotential,slopeDiffs,dists] =...
                findPotentialMergeables(lines,slopeDiff,distDiff,overlapTH)
            %Finds the indices of the lines that are potentially mergeable
            %with the first line of the input array. The first line itself
            %is not included.
            %
            %   The function returns a logical vector defining which lines
            %   are potentially mergeable lines. In addition, the function
            %   returns the slope differences and distances between the
            %   first line and all other lines.
            
            % Find all lines that have a similar slope to the
            % current line
            [locsSlope,slopeDiffs] = LineSearcher.findParallelLines(lines,...
                slopeDiff);
            % Find all lines whose start or end point is close to
            % the start or end point of the current line
            [locsDist,dists] = LineSearcher.findCloseLines(lines,distDiff);
            
            
            % If the current line is a vertical line (abs(slope) > 45),
            % find all lines whose range in the direction of the y-axis
            % does not largely overlap with the y-range of the current line
            if abs(lines(1,6)) > 45
                
                % Get the y-ranges of all lines
                ranges = sort(lines(:,3:4),2);
                
            % If the current line is a horizontal line (abs(slope) <= 45),
            % find all lines whose range in the direction of the x-axis
            % does not largely overlap with the x-range of the current line
            else
                
                % Get the x-ranges of all lines
                ranges = sort(lines(:,1:2),2);
                
            end
            
            locsOverlap = LineSearcher.findNotOverlappingLines(ranges,...
                overlapTH);

            % Find the lines that can potentially be merged with
            % the current line
            locsPotential = locsSlope & locsDist & locsOverlap;
            locsPotential(1) = 0;
        end
        
        function locsOverlap = findNotOverlappingLines(ranges,overlapTH)
            %findNotOverlappingLines Finds the ranges that do not largely
            %overlap with the first range in the given array.
            %   locsOverlap = findNotOverlappingLines(ranges,overlapTH)
            %   Takes two input arguments:
            %   ranges: An nx2 array of ranges. The first row of the array
            %   contains the range to which the other ranges will be
            %   compared.
            %   overlapTH: The threshold overlap percentage determining
            %   whether an overlap percentage is large or not. 
            
            % Go through the range of each line and store the lengths
            % of the "overlap" of each line segment
            overlaps = zeros(size(ranges,1),1);
            for i = 1:size(ranges,1)
                % Partially below the range
                if ranges(i,1) < ranges(1,1) && ranges(i,2) > ranges(1,1)...
                        && ranges(i,2) <= ranges(1,2)
                    overlaps(i) = (ranges(i,2)-ranges(1,1))/min(range(...
                        ranges(i,1:2)),range(ranges(1,1:2)));
                % Fully within the range, but shorter line
                elseif ranges(i,1) >= ranges(1,1) && ranges(i,2) <= ranges(1,2)
                    overlaps(i) = 1;
                % Fully within the range, but longer line
                elseif ranges(i,1) <= ranges(1,1) && ranges(i,2) >= ranges(1,2)
                    overlaps(i) = 1;
                % Partially above the range
                elseif ranges(i,1) >= ranges(1,1) && ranges(i,1)...
                        < ranges(1,2) && ranges(i,2) > ranges(1,2)
                    overlaps(i) = (ranges(1,2)-ranges(i,1))/min(range(...
                        ranges(i,1:2)),range(ranges(1,1:2)));
                end
                % The cases fully below and fully above the range are
                % omitted, as they would result in overlap 0, which is
                % already the default value.
            end
            
            % Get the locations of the ranges with no significant overlap
            % with the first range.
            locsOverlap = overlaps*100 < overlapTH;
        end

        function [locsParallel,slopeDiffs] = findParallelLines(lines,...
                slopeDiff)
            %Finds all the lines whose slope difference between the first
            %line in the input array lines is smaller or equal to the given
            %threshold.
            
            slopeDiffs = min(abs(lines(:,6) - lines(1,6)),...
                (90 - abs(lines(:,6)) ) + (90 - abs(lines(1,6))));
            
            locsParallel = slopeDiffs <= slopeDiff;

        end

        function [locsClose,minDist] = findCloseLines(lines,distDiff)
            %Finds all the lines whose start or end point's distance to the
            %start or end point of the first line in the lines array is
            %smaller or equal to the given threshold.

            % Start to start
            distSS = (lines(:,1)-lines(1,1)).^2 +...
                (lines(:,3)-lines(1,3)).^2;

            % Start to end
            distSE = (lines(:,2)-lines(1,1)).^2 +...
                (lines(:,4)-lines(1,3)).^2;

            % End to start
            distES = (lines(:,1)-lines(1,2)).^2 +...
                (lines(:,3)-lines(1,4)).^2;

            % End to end
            distEE = (lines(:,2)-lines(1,2)).^2 +...
                (lines(:,4)-lines(1,4)).^2;

            % Minimum distance
            minDist = min([distSS,distSE,distES,distEE],[],2);

            % Close points
            locsClose = minDist <= distDiff;
        end
        
        function lines = mergeAndRemove(lines,I)
            % Merges the first line of the input array with the line at
            % position I (I is a logical index vector, not a linear index).
            % Removes the original lines from the input array and stores
            % the merged line as the first line of the input array.
            
            % Find the coordinate combination that results in the longest
            % line segment
            newLines = [lines(1,1), lines(I,2),lines(1,3), lines(I,4);...
                        lines(1,1), lines(I,1),lines(1,3), lines(I,3);...
                        lines(1,2), lines(I,2),lines(1,4), lines(I,4);...
                        lines(1,2), lines(I,1),lines(1,4), lines(I,3)];
            lengths = LineSearcher.lineLength(newLines);
            
            newLine = newLines(lengths == max(lengths),:);
            
            % Make sure the start x coordinate is smaller than the end x
            % coordinate
            newLine = sortrows([newLine(1) newLine(3);newLine(2) newLine(4)]);
            newLine = [newLine(1,1) newLine(2,1) newLine(1,2) newLine(2,2)];
                        
            % Calculate the length and slope of the merged line
            newLine(5) = LineSearcher.lineLength(newLine);
            newLine(6) = LineSearcher.lineSlope(newLine);

            % Remove the original lines from the lines array
            I(1) = 1;
            lines(I,:) = [];

            % Store the merged line as the first line of the lines
            % array
            lines = [newLine;lines];
        end

    end
end













