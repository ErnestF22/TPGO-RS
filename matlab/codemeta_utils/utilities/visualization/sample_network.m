N=7;       % number of nodes

%generate ground truth
Gitruth=zeros(4,4,N);
for(inode=1:N)
    Gitruth(1:3,1:3,inode)=abspose2rot(0,pi/8*mod(inode-1,2)+pi/2,(inode-1)/N*2*pi);
    Gitruth(1:3,4,inode)=Gitruth(1:3,1:3,inode)*[0;0;-8];
    Gitruth(4,4,inode)=1;
end

%generate structure and images
sigmanoisePoints=1; %in pixels

X=9*(rand(3,30)-0.5);
impixel=1000;
for(inode=1:N)
    ximage(:,:,inode)=persp_project(dehom(invg(Gitruth(:,:,inode))*[X; ones(1,size(X,2))]));
end

figure(1)
for(inode=1:N)
    draw3dcamera(Gitruth(1:3,1:3,inode),Gitruth(1:3,4,inode))
    hold on
    if(inode==1)
        ti=Gitruth(1:3,4,1);
        text(ti(1),ti(2),ti(3)+1.1,'Camera 1','HorizontalAlignment','Right');
    end
end
plot3(X(1,:),X(2,:),X(3,:),'k.')
hold off

axis equal
view(3)

figure(2)
for(inode=1:N)
    subplot(3,3,inode)
    plot(ximage(1,:,inode),ximage(2,:,inode),'.')
end
%set(gcf,'Renderer','Painters')
