function [nn,ca] =...
    trainConnectedComponentClassifier(DWStats,nonDWStats)
%trainConnectedComponentClassifier Trains connected component classifiers
%using the given training data.
%   [nn,ca] = trainConnectedComponentClassifier(DWStats,nonDWStats) Trains
%   a feedforward neural network for classifying connected components as
%   either belonging or not belonging to fallen trees.
%
%   Input arguments:
%   DWStats: A structure array containing properties of connected
%   components representing fallen trees. Created with the function
%   collectConnectedComponentTrainingData.
%   nonDWStats: A structure array containing properties of connected
%   components not representing fallen trees. Create with the function
%   collectConnectedComponentTrainingData.
%
%   Output arguments:
%   nn: A trained feedforward neural network
%   ca: The classification accuracy of the neural network,

% Transform data into a format that can be used for training. X contains
% the features and Y contains the corresponding labels.
[X,Y] = transformDataToArrays(DWStats,nonDWStats);

% Split the data into training and validation sets
[Xtrain,Ytrain,Xval,Yval] = splitToTrainAndValSets(X,Y,0.7,0.3);

% Train a feedforward neural network for classifying the two types of
% connected components. Evaluate its performance.
[nn,ca] = trainNN(Xtrain,Ytrain,Xval,Yval);

end

function [X,Y] = transformDataToArrays(DWStats,nonDWStats)
% Transforms the labelled connected components into an array format that
% can be used for training. Extract the relevant features.

% Initialize a vector for storing the label of each connected component (1
% for fallen trees, 0 for non fallen trees)
Y = [ones(length(DWStats),1); zeros(length(nonDWStats),1)];

% Initialize arrays for storing the values of different properties of the
% connected components representing fallen trees
AreaDW = zeros(1,length(DWStats));
MajorAxisLengthDW = zeros(1,length(DWStats));
MinorAxisLengthDW = zeros(1,length(DWStats));
EccentricityDW = zeros(1,length(DWStats));
ConvexAreaDW = zeros(1,length(DWStats));
CircularityDW = zeros(1,length(DWStats));
FilledAreaDW = zeros(1,length(DWStats));
EulerNumberDW = zeros(1,length(DWStats));
EquivDiameterDW = zeros(1,length(DWStats));
SolidityDW = zeros(1,length(DWStats));
ExtentDW = zeros(1,length(DWStats));
PerimeterDW = zeros(1,length(DWStats));
PerimeterOldDW = zeros(1,length(DWStats));
MaxFeretDiameterDW = zeros(1,length(DWStats));
MaxFeretAngleDW = zeros(1,length(DWStats));
MinFeretDiameterDW = zeros(1,length(DWStats));
MinFeretAngleDW = zeros(1,length(DWStats));

% Store the properties of the connected components into the arrays
for i = 1:length(DWStats)

    AreaDW(i) = DWStats(i).Area;
    MajorAxisLengthDW(i) = DWStats(i).MajorAxisLength;
    MinorAxisLengthDW(i) = DWStats(i).MinorAxisLength;
    EccentricityDW(i) = DWStats(i).Eccentricity;
    ConvexAreaDW(i) = DWStats(i).ConvexArea;
    CircularityDW(i) = DWStats(i).Circularity;
    FilledAreaDW(i) = DWStats(i).FilledArea;
    EulerNumberDW(i) = DWStats(i).EulerNumber;
    EquivDiameterDW(i) = DWStats(i).EquivDiameter;
    SolidityDW(i) = DWStats(i).Solidity;
    ExtentDW(i) = DWStats(i).Extent;
    PerimeterDW(i) = DWStats(i).Perimeter;
    PerimeterOldDW(i) = DWStats(i).PerimeterOld;
    MaxFeretDiameterDW(i) = DWStats(i).MaxFeretDiameter;
    MaxFeretAngleDW(i) = DWStats(i).MaxFeretAngle;
    MinFeretDiameterDW(i) = DWStats(i).MinFeretDiameter;
    MinFeretAngleDW(i) = DWStats(i).MinFeretAngle;
    
end

% Initialize arrays for storing the values of different properties of the
% connected components not representing fallen trees
AreaNonDW = zeros(1,length(nonDWStats));
MajorAxisLengthNonDW = zeros(1,length(nonDWStats),1);
MinorAxisLengthNonDW = zeros(1,length(nonDWStats),1);
EccentricityNonDW = zeros(1,length(nonDWStats));
ConvexAreaNonDW = zeros(1,length(nonDWStats));
CircularityNonDW = zeros(1,length(nonDWStats));
FilledAreaNonDW = zeros(1,length(nonDWStats));
EulerNumberNonDW = zeros(1,length(nonDWStats));
EquivDiameterNonDW = zeros(1,length(nonDWStats));
SolidityNonDW = zeros(1,length(nonDWStats));
ExtentNonDW = zeros(1,length(nonDWStats));
PerimeterNonDW = zeros(1,length(nonDWStats));
PerimeterOldNonDW = zeros(1,length(nonDWStats));
MaxFeretDiameterNonDW = zeros(1,length(nonDWStats));
MaxFeretAngleNonDW = zeros(1,length(nonDWStats));
MinFeretDiameterNonDW = zeros(1,length(nonDWStats));
MinFeretAngleNonDW = zeros(1,length(nonDWStats));

% Store the properties of the connected components into the arrays
for i = 1:length(nonDWStats)

    AreaNonDW(i) = nonDWStats(i).Area;
    MajorAxisLengthNonDW(i) = nonDWStats(i).MajorAxisLength;
    MinorAxisLengthNonDW(i) = nonDWStats(i).MinorAxisLength;
    EccentricityNonDW(i) = nonDWStats(i).Eccentricity;
    ConvexAreaNonDW(i) = nonDWStats(i).ConvexArea;
    CircularityNonDW(i) = nonDWStats(i).Circularity;
    FilledAreaNonDW(i) = nonDWStats(i).FilledArea;
    EulerNumberNonDW(i) = nonDWStats(i).EulerNumber;
    EquivDiameterNonDW(i) = nonDWStats(i).EquivDiameter;
    SolidityNonDW(i) = nonDWStats(i).Solidity;
    ExtentNonDW(i) = nonDWStats(i).Extent;
    PerimeterNonDW(i) = nonDWStats(i).Perimeter;
    PerimeterOldNonDW(i) = nonDWStats(i).PerimeterOld;
    MaxFeretDiameterNonDW(i) = nonDWStats(i).MaxFeretDiameter;
    MaxFeretAngleNonDW(i) = nonDWStats(i).MaxFeretAngle;
    MinFeretDiameterNonDW(i) = nonDWStats(i).MinFeretDiameter;
    MinFeretAngleNonDW(i) = nonDWStats(i).MinFeretAngle;
    
end

X = [AreaDW AreaNonDW;...
    CircularityDW CircularityNonDW;...
    ConvexAreaDW./AreaDW ConvexAreaNonDW./ConvexAreaNonDW;...
    EccentricityDW EccentricityNonDW;...
    ExtentDW ExtentNonDW;...
    FilledAreaDW./AreaDW FilledAreaNonDW./AreaNonDW;...
    MajorAxisLengthDW./AreaDW MajorAxisLengthNonDW./AreaNonDW;...
    MaxFeretDiameterDW./AreaDW MaxFeretDiameterNonDW./AreaNonDW;...
    MinorAxisLengthDW./AreaDW MinorAxisLengthNonDW./AreaNonDW;...
    PerimeterDW./AreaDW PerimeterNonDW./AreaNonDW;...
    SolidityDW SolidityNonDW
    MinFeretDiameterDW./AreaDW MinFeretDiameterNonDW./AreaNonDW]';

end

function [Xtrain,Ytrain,Xval,Yval] = splitToTrainAndValSets(X,Y,tProp,...
    vProp)

% Split the dataset into training and validation sets
[IdxTrain,IdxVal,~] = splitToTrainValTest(randperm(size(X,1))',tProp,...
    vProp,0);
Xtrain = X(IdxTrain,:);
Ytrain = Y(IdxTrain);
Xval = X(IdxVal,:);
Yval = Y(IdxVal);

end

function [nn,CA] = trainNN(Xtrain,Ytrain,Xval,Yval)
% Trains a neural network and calculates its classification accuracy on
% validation data.

net = patternnet();
net.divideParam.trainRatio = 1;
net.divideParam.valRatio = 0;
net.divideParam.testRatio = 0;
nn = train(net,Xtrain',Ytrain');

% Classification accuracy
predObs = [(round(nn(Xval')))',Yval];
CA = classificationAccuracy(predObs);
end

function CA = classificationAccuracy(predObs)
% Calculates the classification accuracy from an array in which the first
% column contains the predictions and the second column contains the
% corresponding observations.

CA = sum(predObs(:,1) == predObs(:,2))/size(predObs,1);
end

