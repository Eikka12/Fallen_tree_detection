function gridsArray = processSegments(pointsCell,cellSize,targetSize)
    % processSegments Processes the given segments to the input format of a
    % convolutional neural network.
    %   processSegments(pointsCell,cellSize,targetSize) Takes three
    %   input arguments:
    %   pointsCell: A cell array containing the laser points of each
    %   segment to be processed.
    %   cellSize: The cell size of the grids created from each segment.
    %   targetSize: A 2-dimensional vector specifying the number of rows
    %   and columns in the output grids.
    %
    %   The function processes the input segments to the input format of a
    %   convolutional neural network. The processing steps are:
    %   1. Create a grid from each segment.
    %   2. Scale the grid so that either the number of rows or columns is
    %   equal to the number of rows or columns defined by targetSize.
    %   3. Pad the array with zeros so that the resulting grid is of size
    %   targetSize(1)xtargetSize(2).
    %   4. Create an targetSize(1)xtargetSize(2)x3 grid by copying
    %   the grid twice. This step is performed, as the convolutional neural
    %   networks requires three channels (RGB) as input.
    %   5. Add all processed grids together to form a 4-dimensional array
    %   of size targetSize(1)xtargetSize(2)x3xN, where N is the
    %   number of grids.
    
    % Create grids from each segment
    gridsCell = createGrids(pointsCell,cellSize);
    
    % Create a grid array to which the processed grids will be stored
    gridsArray = zeros(targetSize(1),targetSize(2),3,length(...
        gridsCell));
    
    % Process each grid
    for i = 1:length(gridsCell)
        gridsArray(:,:,:,i) = processGrid(gridsCell{i},targetSize);
    end
end

function gridsCell = createGrids(segments,cellSize)
    % Creates a rasterized representation of each segment and returns the
    % representations as a cell array in which each cell contains the grid
    % of one segment.
    
    % Create a cell array for storing the grids
    gridsCell = cell(size(segments));
    
    % Go through each segment
    for i = 1:length(segments)
        % Create a grid representation of the segment and store it
        raster = Grid(segments{i},cellSize);
        raster.cellCount(1);
        [gridsCell{i},~] = raster.getGrid;
    end
end

function processedGrid = processGrid(g,targetSize)
    % processGrid Processes one grid
    
    % If the grid contains no cells
    if isempty(g)
        % Create a grid with one cell with the value zero
        g = zeros(1,1);
    end
    
    % Set all non-zero cells to 255
    g(g~=0) = 255;
    
    % Get the grid size
    sz = size(g);
   
    % Scale the grid to fit the targetSize as well as possible
    scaleR = targetSize(1)/sz(1);
    scaleC = targetSize(2)/sz(2);

    if scaleR <= scaleC
        g = imresize(g,[targetSize(1) NaN]);
    else
        g = imresize(g,[NaN targetSize(2)]);
    end
    
    % Get the new grid dimensions
    sz = size(g);

    % Pad the grid so that it becomes the same size as the target size
    g = padarray(g,[targetSize(1)-sz(1),...
        targetSize(2)-sz(2)],'post');
    
    % Copy the grid to create a three-dimensional grid
    processedGrid(:,:,1) = g;
    processedGrid(:,:,2) = g;
    processedGrid(:,:,3) = g;
    
    
end
    