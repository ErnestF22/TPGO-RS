function notes_plots
L=10;
Nw=6;

%w=cnormalize(randn(2,1));
%b=2*L*(rand()-0.5);

wVert=[1;0];
wHoriz=[0; 1];

w=[repmat(wVert,1,Nw/2) repmat(wHoriz,1,Nw/2)];
v=repmat(linspace(-L,L,Nw/2),1,2);

% w=[
%     1  0;
%     0 -1;
%   ];
% v=[0 0];

%a=[0;0;1;1];


NPoints=30;
%xPositive=rand(2,NPoints);
%xNegative=[-rand(1,2*NPoints) rand(1,NPoints); rand(1,NPoints) -rand(1,2*NPoints)];
xxBase=linspace(-0.2,0.2,round(sqrt(NPoints)))
[xBaseX,xBaseY]=meshgrid(xxBase,xxBase);
xBaseX=reshape(xBaseX,1,[]);
xBaseY=reshape(xBaseY,1,[]);
t=0.3;
xPositive=[xBaseX;xBaseY]+t;
xNegative=[ xBaseX+t xBaseX-t xBaseX-t;
            xBaseY-t xBaseY-t xBaseY+t];

x=[xPositive xNegative];
y=[+ones(1,size(xPositive,2)) -ones(1,size(xNegative,2))]';

% figure(1)
% showHingeResponse(w(:,1),b(1))

[a,b]=hingeSparseBoost(x,y,w,v,'C',100);

figure(2)
showHingeEnsambleResponse(w,v,a,b,'displayContour')
hold on
plotPoints(xPositive,'r*')
plotPoints(xNegative,'b*')
hold off

figure(3)
subplot(2,1,1)
plot(hingeEnsambleResponse(x,w,v,a,b))
subplot(2,1,2)
plot(a)

disp([w -w;v -v;a'])
