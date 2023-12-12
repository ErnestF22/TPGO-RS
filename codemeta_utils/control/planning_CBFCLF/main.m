%%function [K1]= main
 close all
 clear 
[n,Ah,bh,y,Cl,Cb,Wl,Wb,Au,bu,Ax,bx]=data;
yt=y;
plot_polygon(y,'K')%plots the polygon
plot(yt(1,:),yt(2,:),'r*')%plots the landmarks
%finding K1:

 [K1] =  optFirstorder(Wb,Wl,Cb,Cl,n,Ah,Au,Ax,bx,bh,bu,yt);
%[K1] =  optFirstorderNOTvectorized(Wb,Wl,Cb,Cl,n,Ah,Au,Ax,bx,bh,bu,yt);

plotController(yt,y,K1,'b')%plots the control flow
% x0=yt(:,1)+[-10;15];
% simulationOFpath(K1,yt,x0,0.1,'r*')%plots the path
% end



%%
return
yt=transferedLandmarks(y,30);
% yt=10*y+[25;300]*[1 1 1 1];
[result2] = optFirstorderNOTvectorized(Wb,Wl,Cb,Cl,n,Ah,Au,Ax,bx,bh,bu,yt);
K1=result2.K;
plotController(yt,y,K1,'b')%plots the control flow
plot_polygon(y,'K')%plots the polygon
plot(yt(1,:),yt(2,:),'r*')%plots the landmarks


k1=result1.K
(y-[20;60]*[1 1 1 1])*k1
sum(k1)
k2=result2.K
(yt-[20;60]*[1 1 1 1])*k2
sum(k2)

% a=R*y;
% a=[a;1 1 1 1];
% k22=(a'*a)^-1*a'*[y;1 1 1 1]*k1;