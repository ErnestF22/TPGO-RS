function subData=ransacFunSelectionColumns(data,idx)
if nargin==1
    subData=size(data,2);
else
    subData=data(:,idx,:);
end
