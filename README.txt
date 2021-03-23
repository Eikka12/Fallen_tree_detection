This set of folders contains the source code for the ALS-based fallen tree 
detection algorithm developed by Heinaro et al. 2021. The algorithm is 
used via the function findFallenTrees that is located in the main folder. 
A short example usage of the function can be found from the example folder.
Check the function documentation for a more detailed description of the 
inputs, outputs and usage of the function. NOTE: The script startup.m must
be called before running the function, as it adds all the required 
filepaths to the MATLAB path.

The fallen tree detection algorithm consists of four steps:

1. Reading and preprocessing the data
2. Filtering the point cloud using classification based on connected 
component analysis
3. Detecting fallen trees using iterative Hough transform -based line 
detection
4. Removing false fallen tree segments using a convolutional neural network

Steps 2 and 4 are optional. By default, the function findFallenTrees
performs step 2 using a default classifier that can be found from the 
folder trained_classifiers, whereas step 4 is not performed. Steps 2 and 4 
can use user-defined classifiers trained using the functions found in the 
folders connected_component_training (step 2) and final_classifier_training
(step 4).
