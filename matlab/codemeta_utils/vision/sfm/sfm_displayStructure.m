function sfm_displayStructure(data)
if isfield(data,'pose')
    poseMemberName='pose';
elseif isfield(data,'poseTruth')
    poseMemberName='poseTruth';
    display('Displaying ground truth pose')
else
    error('Could not find a pose member')
end
testNetworkDisplay(data.(poseMemberName),'references')
hold on
plotPoints(data.structureFiltered.location);
hold off
axis square
axis equal
grid on
