%Return the result of now as a filename-friendly string
function s=nowFileName()
s=regexprep(datestr(now),'[ :]','_');
