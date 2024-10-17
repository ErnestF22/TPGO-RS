function testNetworkShowErrors(t_node,varargin)
flagHasTruth=isfield(t_node,'gitruth');
methodAbsolutePoses='reference';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

figure(11)
set(gcf,'Name','Network')
if flagHasTruth
    t_node=testNetworkCompensate(t_node,'Poses');
    testNetworkDisplay(t_node,'DisplayEdges','methodabsoluteposes',methodAbsolutePoses);
    hold on
end
testNetworkDisplay(t_node,'Member','gi','Estimated','methodabsoluteposes',methodAbsolutePoses);
if flagHasTruth
    hold off
end

[rotErr,translErr,scaleRatios]=testNetworkComputeErrors(t_node,'methodabsoluteposes',methodAbsolutePoses);
if flagHasTruth
    [rotErrMs,translErrMs,scaleRatiosMs]=testNetworkComputeErrors(t_node,'Measurements','methodabsoluteposes',methodAbsolutePoses);
end

%Note:
% - measurements errors should be all zero
% - errors should be in the order of 2 times the variances given during
%   state initialization
figure(12)
set(gcf,'Name','Rotation and Translation errors')
switch flagHasTruth
    case true
        subplot(2,1,1)
        hist([rotErr rotErrMs]*180/pi)
        legend('Estimated VS Ground Truth', 'Measurements VS Ground Truth')
        subplot(2,1,2)
        hist([translErr translErrMs]*180/pi)
        disp('Estimated VS Ground Truth')
        displayStatRTError(rotErr,translErr)
        disp('Measurements VS Ground Truth')
        displayStatRTError(rotErrMs,translErrMs)
        disp(['Geometric std of ratios (est/mes)', num2str(geostd(scaleRatios)), '/', num2str(geostd(scaleRatiosMs))])
    case false
        subplot(2,1,1)
        hist(rotErr*180/pi)
        legend('Estimated VS Measurements')
        subplot(2,1,2)
        hist(translErr*180/pi)
        disp('Estimated VS Measurements')
        displayStatRTError(rotErr,translErr)
end
subplot(2,1,1)
title('Rotation angle Errors')
xlabel('Degrees')
subplot(2,1,2)
title('Translation Angle Errors')
xlabel('Degrees')


function displayStatRTError(rotErr, translErr)
%NOTE: errors must be passed in radians but are shown in degrees
fprintf('\tRot. angle errors [deg]\n')
displayStatError(rotErr*180/pi)
fprintf('\tTransl. angle errors [deg]\n')
displayStatError(translErr*180/pi)

function displayStatError(err)
fprintf('\t\tMean  :%.4f\n', mean(err))
fprintf('\t\tStd   :%.4f\n', std(err))
fprintf('\t\tMedian:%.4f\n', median(err))
