function plotPointsNumbers(xPoints,varargin)
for iPoint=1:size(xPoints,2)
    text(xPoints(1,iPoint),xPoints(2,iPoint),num2str(iPoint),varargin{:})
end