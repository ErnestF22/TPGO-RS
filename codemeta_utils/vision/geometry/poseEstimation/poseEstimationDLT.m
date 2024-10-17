%function [Rest,Test]=pose_estimation(x_model,p)
%Estimate the pose of an object having 3D shape x_model from the image
%points p
function GEst=poseEstimationDLT(x_model,p,varargin)

flagnormalize=true;
flagsecondorder=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'nonormalize'
            flagnormalize=false;
        case 'secondorder'
            flagsecondorder=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

p=homogeneous(p,3);
n=size(x_model,2);    %number of points in the model

if(size(p,1)==2)
    p=[p;ones(1,n)];
end
if(size(x_model,1)==3) 
    x_model=[x_model;ones(1,n)];
end

if(flagnormalize)
    [p,Tnorm_p]=normalize_anisot(p);
    [x_model,Tnorm_X]=normalize_anisot(x_model);
end

%find projection matrix P using hat(p)Px_model=0
A=[];
for(i=1:n)
    Ai=kron(x_model(:,i)',hat(p(:,i)));
    A=[A;Ai(1:2,:)];
end

if(flagsecondorder)
    for i=1:n
        for j=1:n
            if i~=j
                Ai1=kron((x_model(:,i))',(hat(p(:,i))*p(:,j))');
                Ai2=kron((x_model(:,i))',(hat(p(:,i))*p(:,j))');
                A=[A;Ai1;Ai2];
            end
        end
    end
end    

[U,S,V]=svd(A,'econ');

P=reshape(V(:,end),3,4);

if(flagnormalize)
    P=inv(Tnorm_p)*P*Tnorm_X;
end

Xc=P*mean(x_model,2);
if(Xc(3)<0)
    P=-P;
end

%decompose P in rotation and translation
%compute RQ decomposition from the SVD of P(:,1:3)
[Up,Sp,Vp]=svd(P(:,1:3));

if(det(P(:,1:3))>0)
    REst=Up*Vp';
else
    REst=Up*[1 0 0;0 1 0; 0 0 -1]*Vp';
end
K=P(:,1:3)*REst';

scale=K(3,3);
K=K/scale;

TEst=P(:,4)/scale;

GEst=RT2G(REst,TEst);

%normalize homogeneous points
function [p1,T]=normalize_anisot(p)
D=size(p,1);
T=eye(D);
for(d=1:D-1)
    m=mean(p(d,:));
    v=var(p(d,:)-m);
    T(d,d)=1/sqrt(v);
    T(d,D)=-m;
end
p1=T*p;



