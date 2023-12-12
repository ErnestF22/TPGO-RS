%Evaluate empirical cumulative distribution
%function [d,x]=cumDistEval(data,x)
%Returns value of empirical distribution d of the data evaluated at
%locations x. More in detail, d(i) contains the fraction (between 0 and 1)
%of values in data less than or equal to x(i)
function [d,x]=cumDistEval(data,x)
NData=length(data);
c=cumCountEval(data,x);
d=c/NData;


