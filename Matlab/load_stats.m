function data=load_stats(fileName)

AllData = importdata(fileName);

% Create new variables in the base workspace from those fields.
%vars = fieldnames(AllData);

data=AllData.data;

return