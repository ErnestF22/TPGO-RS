function fontSizeSetForPublication(fontSize)
warning('Deprecated. Use setFigFontSize before opening the figure instead')
if ~exist('fontSize','var') || isempty(fontSize)
    fontSize=8;
end
set(get(gca,'xlabel'),'FontSize',fontSize)
set(get(gca,'ylabel'),'FontSize',fontSize)
set(gca,'FontSize',fontSize)
set(legend,'FontSize',fontSize)
