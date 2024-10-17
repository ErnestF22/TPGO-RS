%function multiCalibration_posesDisplay(posesAbsolute)
%Display absolute poses, and use different shapes according to
%poseAbsolute.nodeName
function multiCalibration_posesDisplay(posesAbsolute)
methodAbsolutePoses='reference';
NNodes=length(posesAbsolute);
for iNode=1:NNodes
    G=posesAbsolute(iNode).G;
    switch posesAbsolute(iNode).nodeName
        case 'velodyne'
            draw3dcameraFromG(G,'shape','velodyne','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
        case 'hokuyo'
            draw3dcameraFromG(G,'shape','hokuyo','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
        otherwise
            draw3dcameraFromG(G,'scale',0.25,'methodAbsolutePoses',methodAbsolutePoses);
    end
    hold on
end
hold off
axis equal
