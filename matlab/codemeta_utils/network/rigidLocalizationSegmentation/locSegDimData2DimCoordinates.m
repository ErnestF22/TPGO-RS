function [dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData)
switch dimData
    case 2
        dimCoordinates=3;
        dimRotations=1;
    case 3
        dimCoordinates=6;
        dimRotations=3;
end
