%Returns true if it is a plot style string
%function flag=isStyleString(s)
function flag=isStyleString(s)
flag=~isempty(regexpi(s,'[bgrcmykwoxsdvph\.\-:<>\^\*\+]*$'));
