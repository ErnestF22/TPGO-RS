function status = rigidDyn_odeplot(t,y,flag,varargin)
s=zeros(4,1);
%[R,w,T,v]=rigidDyn_stateUnpackRT(y);
figure(1)
s(1)=odeplot(t,y,flag,varargin{:});
figure(2)
s(2)=odeplot(t,y,flag,varargin{:});
figure(3)
s(3)=odeplot(t,y,flag,varargin{:});
figure(4)
s(4)=odeplot(t,y,flag,varargin{:});
status=max(s);


