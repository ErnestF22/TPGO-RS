function [geoTagIdx,geoTag]=testNetworkGetGeoTags(t_node)

structType=testNetworkDetectStructType(t_node);

if ~isfield(t_node,'flagHasGeoTag')
    geoTagIdx=[];
    geoTag=[];
else
    switch structType
        case 'single'
            geoTagIdx=find(t_node.flagHasGeoTag);
            geoTag=t_node.geoTag(:,geoTagIdx);
        case 'array'
            geoTagIdx=find([t_node.flagHasGeoTag]);
            geoTag=cat(2,t_node(geoTagIdx).geoTag);
    end
    geoTag=geoTag(1:2,:);
end
