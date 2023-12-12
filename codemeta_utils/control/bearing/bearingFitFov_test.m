function bearingFitFov_test
close all

allAngles=[pi/4 3*pi/4];
allFovs=1.8*allAngles;
N=30;
z=zeros(1,N);

cnt=1;
for iAngle=1:2
    for iFov=1:2
        a=allAngles(iAngle);
        fov=allFovs(iFov);
        
        t=2*a*linspace(0,1,N).^2;
        y=[cos(t);sin(t)];
        yHeading=bearingFitFov(y);
        subplot(2,2,cnt)
        bearingPlotFov([0;0],yHeading,fov)
        hold on
        quiver(z,z,y(1,:),y(2,:),1);
        hold off
        axis equal
        
        cnt=cnt+1;
    end
end