function [dTrain,dVal,dTest] = splitToTrainValTest(d,tr,v,te)
%splitToTrainValTest Splits the given data into training, validation and
%test data.
%   [dTrain,dVal,dTest] = splitToTrainValTest(d,tr,v,te) Takes four input
%   arguments:
%   d: The data to be splitted. The rows of the data must contain the
%   individual observations.
%   tr: The number or fraction of training data observations to be
%   extracted
%   v: The number or fraction of validation data observations to be
%   extracted
%   te: The number or fraction of test data observations to be extracted
%
%   The function splits the given data into training, validation and test
%   sets. If tr, v and te are values larger than one, they are interpreted
%   as numbers of observations in each set. The sum of the numbers of
%   observations must not exceed the number of observations in the data. If
%   tr, v and te are values smaller than one, they are interpreted as
%   fractions of the data. The sum of the fractions must not exceed 1.

% If all values are smaller than or equal to one and the sum of the values
% is smaller than or equal to 1
if sum([tr,v,te]) <= 1
    % Interpret the values as fractions and calculate the corresponding
    % numbers of observations
    tr = round(tr*size(d,1));
    v = floor(v*size(d,1));
    te = floor(te*size(d,1));
    
% If the input is not valid 
elseif sum([tr,v,te]) > size(d,1)
    error('Input is not valid. Give tr, v and te as fractions or numbers and check that the sum of tr, v and te does not exceed the maximum value.')
end

% Select the training data
[dTrain,dRemaining] = selectNObservations(d,tr);

% Select the validation data
[dVal,dRemaining]= selectNObservations(dRemaining,v);

% Select the test data
[dTest,~] = selectNObservations(dRemaining,te);

end

function [extracted,remaining] = selectNObservations(d,N)
%selectNObservations Selects N observations from the data randomly.
%   [extracted,remaining] = selectNObservations(d,N) Takes two input
%   arguments:
%   d: A dataset containing observations at each row.
%   N: The number of observations to be extracted
%
%   The function extracts N random observations from the data. It then
%   returns the extracted observations as well as the remaining
%   observations as two separate arrays.

% Get the indices of the observations to be extracted
idx = randperm(size(d,1),N);
% Convert the indices into a logical vector to enable inverse queries
locs = false(size(d,1),1);
locs(idx) = 1;

% Extract the observations
extracted = d(locs);
% Extract the remaining observations
remaining = d(~locs);

end