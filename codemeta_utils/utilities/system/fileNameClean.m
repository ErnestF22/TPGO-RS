%Replaces the characters ' :.' with _
%function s=fileNameClean(s)
function s=fileNameClean(s)
s=regexprep(s,'[ \.:]','_');