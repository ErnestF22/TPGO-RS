%function prevSize=setFigFontSize(newSize)
%Sets the default font size for plots to NEWSIZE pt. Returns the current
%font size
function prevSize=setFigFontSize(newSize)

if nargout>0
    prevSize=get(0, 'DefaultAxesFontSize');
end
set(0, 'DefaultAxesFontSize', newSize)


