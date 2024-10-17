%function prevFont=setFigFont(newFont)
%Sets the default font for plots. Returns the current font
function prevFont=setFigFont(newFont)

if nargout>0
    prevFont=get(0, 'DefaultAxesFontName');
end
set(0, 'DefaultAxesFontName', newFont)


