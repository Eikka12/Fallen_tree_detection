This folder contains two functions:

- collectConnectedComponentTrainingData
- trainConnectedComponentClassifier

These two functions are part of the cleaning algorithm that removes groups of points that do not represent fallen trees from a laser point cloud. The first function is used for collecting training data for connected component classification. The second function uses this training data to train a classifier that classifies connected components as either representing or not representing parts of fallen trees. See the documentation of the functions for a more detailed description of the inputs and outputs as well as how the functions work. 