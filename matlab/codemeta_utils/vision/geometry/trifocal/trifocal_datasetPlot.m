function trifocal_datasetPlot
data=trifocal_dataset();
testNetworkDisplay(data.G)
hold on
NLines=size(data.l3d,3);
for iLine=1:NLines
    draw3dLine(data.l3d(:,:,iLine),'side',3,'plotOptions',{'color',rand(1,3)});
end
hold off
axis tight
