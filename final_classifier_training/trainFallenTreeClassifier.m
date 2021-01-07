function trainedNet = trainFallenTreeClassifier(segments,labels,cellSize)
%trainFallenTreeClassifier Trains a fallen tree classifier
%   trainedNet = trainFallenTreeClassifier(segments,labels,cellSize) Takes
%   three input arguments:
%   segments: A cell array containing point cloud representations of fallen
%   tree segments. The training segments can be extracted using the
%   function collectTrainingSegments.
%   labels: A logical vector containing the labels of the segments. 1 means
%   that the segment represents a fallen tree and 0 means that the segment
%   does not represent a fallen tree. The labels can be extracted using the
%   function collectTrainingSegments.
%   cellSize: The cell size of the binary image created from each segment.
%
%   The function trains a convolutional neural network for classifying
%   fallen tree segments as either real or false fallen tree segments. The
%   classifier is a tuned version of AlexNet - a pretrained convolutional
%   neural network.

%% Equalize the number of training samples in both classes
numNoTree = sum(labels == 0);
numTree = sum(labels == 1);

% If the number of samples in class 0 is larger
if numNoTree > numTree
    noTreeLocs = find(labels == 0);
    notIncluded = noTreeLocs(randperm(numNoTree,numNoTree-numTree));

% If the number of samples in class 1 is larger
elseif numTree > numNoTree
    treeLocs = find(labels == 1);
    notIncluded = treeLocs(randperm(numTree,numTree-numNoTree));
end

% Remove some of the samples belonging to the larger class
segments(notIncluded) = [];
labels(notIncluded) = [];

%% Process the training segments to the input format required by alexnet
% The input format is a 227x227x3xN array, an array containing N 227x227x3
% images, each representing one segment)
images = processSegments(segments,cellSize,[227 227]);

%% Split the images into training and validation images
% Extract the indices of the training data (70 %) and validation data (30
% %)
indices = (1:size(images,4))';
[iTrain,iVal,~] = splitToTrainValTest(indices,0.7,0.3,0);
% Extract the training and validation data
trainingImages = images(:,:,:,iTrain);
validationImages = images(:,:,:,iVal);
trainingLabels = categorical(labels(iTrain));
validationLabels = categorical(labels(iVal));

%% Prepare the pretrained network for training
% Load AlexNet (Requires that the alexnet add-on has been downloaded)
net = alexnet;
layers = net.Layers;

% Replace the final fully connected layer with a new fully connected layer
% with the same size as the number of classes (2). Increase the learning
% rate factors to increase the learning speed of the new layer
newFCLayer = fullyConnectedLayer(2,'Name','new_fc',...
    'WeightLearnRateFactor',10,'BiasLearnRateFactor',10);
layers(23) = newFCLayer;
% Replace the original classification layer with a new classification layer
newClassLayer = classificationLayer('Name','new_classoutput');
layers(25) = newClassLayer;

%% Train the network
% Define the training options
options = trainingOptions('sgdm', ...
    'MiniBatchSize',128, ...
    'MaxEpochs',50, ...
    'InitialLearnRate',1e-4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{validationImages,validationLabels}, ...
    'ValidationFrequency',5, ...
    'ValidationPatience',5, ...
    'Verbose',false, ...
    'Plots','training-progress');

% Train the network
trainedNet = trainNetwork(trainingImages,trainingLabels,layers,options);


end

