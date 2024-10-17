function presentation_sfmExample
load sfmdata_castle
sfm_displayStructure(data)
view(2)
axis([-30 30 -30 20])
set(gca,'Visible','off')
savefigure('~/Documents/repository/JHU/Presentations/Figures/sfmExample','epsc',[400 300]*0.75,2)
