function T=realDyn_integrateAcceleration(t,a,T0,v0)
NIt=length(t);
dt=diff(t);
T=zeros(3,NIt);
T(:,1)=T0;
v=zeros(3,NIt);
v(:,1)=v0;
for it=2:NIt
    v(:,it)=v(:,it-1)+dt(it-1)*a(:,it-1);
    T(:,it)=T(:,it-1)+dt(it-1)*v(:,it-1);
end
