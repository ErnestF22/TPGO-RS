
function [Ropt,tglobalopt]=geod_mean_rot(allR,R0,v)
v=v/norm(v);
% 
% T=linspace(0,2*pi,200);
% d=[];
% for(t=T)
%     d(end+1)=sumdist(t,allR,R0,v);
% end
% plot(T,d)
% hold on

N=size(allR,3);
t=zeros(1,N);
for(r=1:N)
    [S,t(r)]=geodRotProj(allR(:,:,r),R0,v);
end
t=sort(mod(t+pi,2*pi));
% 
% plot(t,zeros(size(t)),'rd');

if(N==1)
    tglobalopt=t;
else
    topt=zeros(1,N);
    topt(1)=mod(fminbnd(@(x) sumdist(x,allR,R0,v),t(N),t(1)+2*pi),2*pi);
    for(r=2:N)
        topt(r)=fminbnd(@(x) sumdist(x,allR,R0,v),t(r-1),t(r));
    end

    dopt=zeros(1,N);
    for(r=1:N)
        dopt(r)=sumdist(topt(r),allR,R0,v);
    end
%     plot(topt,dopt,'ro');
% 
    [mind,minidx]=min(dopt);
    tglobalopt=topt(minidx);
end

Ropt=R0*rot(tglobalopt*v);

function d=sumdist(t,allR,R0,v)
d=sum(rot_dist(allR,R0*rot(t*v)).^2);
    