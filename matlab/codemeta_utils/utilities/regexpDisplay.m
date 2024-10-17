%Display a string together with regular expression matches
function regexpDisplay(str,expression,indent)
if ~exist('indent','var')
    indent='';
end

if isnumeric(expression)
    idx=expression;
else
    idx=regexp(str,expression);
end

if isempty(idx)
    disp('No matches')
else
    indicators=repmat(' ',1,length(str));
    indicators(idx)='^';
    fprintf('%s%s\n%s%s\n',indent,str,indent,indicators)
end
