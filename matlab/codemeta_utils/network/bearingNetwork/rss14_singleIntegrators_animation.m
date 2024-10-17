function rss14_singleIntegrators_animation
close all
fs=setFigFontSize(12);
fn=setFigFont('Verdana');
costNameBearings='cosine';
flagUseRanges=true;
flagRecordVideo=true;
deltaT=10; %step in simulated time between frames
frameRate=30; %video framerate

baseFileName=['bearingNetwork_' costNameBearings '_'];
if ~flagUseRanges
    disp('## Pure bearing formation')
    baseFileName=[baseFileName 'b_'];
else
    disp('## Bearing+distance formation')
    baseFileName=[baseFileName 'bd_'];
end

load([baseFileName 'data'],'t','x0','x','xg','xFinal','t_node',...
    'phi','funsBearings','c','m','d')
if flagUseRanges
    load([baseFileName 'data'],'funsRanges','q')
end

figure(1)
maxT=t(end);
%maxT=120;

if flagRecordVideo
    vWrt=VideoWriter([mfilename '.avi']);
    vWrt.FrameRate=frameRate;
    open(vWrt);
end

figure(1)
set(gcf,'Color',[1 1 1])
for tFrame=1:deltaT:maxT;
    [~,idxT]=min(abs(t-tFrame));
    plot(x0(1,:),x0(2,:),'rx')
    hold on
    plot(squeeze(x(1,:,1:idxT))',squeeze(x(2,:,1:idxT))','-','color',[1 0.75 0])
    plot([xg(1,:); x(1,:,idxT)],[xg(2,:); x(2,:,idxT)],'k:')
    t_node.Ti=x(:,:,idxT);
    bearingNetworkPlot(t_node,'flagPlotRanges',flagUseRanges)
    hold off
    axis([-11 10 -8.5 7.5])
    %title(['t=' num2str(tFrame)])
    pause(0.0001)
    if flagRecordVideo
        writeVideo(vWrt,getframe(gcf));
    end
end

if flagRecordVideo
    close(vWrt);
end

setFigFontSize(fs);
setFigFont(fn);
