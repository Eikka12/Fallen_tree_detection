function class = predictClass(stats,net)
%predictClass Classifies the input connected component into a component
%belonging or not belonging to a fallen tree.
%   class = predictClass(stats) Takes two input arguments: 
%   stats: a structure array that contains statistics computed with
%   regionprops
%   net: a trained neural network used for classification
%
%   The function classifies a connected component into a component
%   belonging (1) or not belonging (0) to a fallen tree based on the given
%   statistics of the connected component

Area = [stats.Area];

X = [[stats.Area]' [stats.Circularity]' [stats.ConvexArea]'./Area'...
    [stats.Eccentricity]' [stats.Extent]' [stats.FilledArea]'./Area'...
    [stats.MajorAxisLength]'./Area' [stats.MaxFeretDiameter]'./Area'...
    [stats.MinorAxisLength]'./Area' [stats.Perimeter]'./Area'...
    [stats.Solidity]' [stats.MinFeretDiameter]'./Area'];

class = net(X');
end

