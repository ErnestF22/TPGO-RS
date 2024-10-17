function subStrings=regexpExtract(str,expression)
[startIndex,endIndex]=regexp(str,expression);
nbStrings=length(startIndex);
subStrings=cell(1,nbStrings);
for iString=1:nbStrings
    subStrings{iString}=str(startIndex(iString):endIndex(iString));
end
