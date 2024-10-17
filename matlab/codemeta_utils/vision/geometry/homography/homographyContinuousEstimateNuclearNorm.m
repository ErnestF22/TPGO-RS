%Factorize continuous homography from multiple views and planes
%function [w,v,n,output]=homographyContinuousEstimateNuclearNorm(H,E)
%Inputs
%   H   [3 x 3 x NPairs] array of continuous-time homographies
%   E   [NPairs x 2] array of pairs [nFrame nPlane], indicating, for each
%       matrix in H, to which frame/plane pair it corresponds
%Outputs
%   w   [3 x NFrames] array of estimated rotational velocities
%   v   [3 x NFrames] array of estimated translational velocities
%   n   [3 x NPlanes] array of estimated scaled plane normals
%   where NFrames=max(E(:,1)) and NPlanes=max(E(:,2))
%
%This function is based on the use of the nuclear norm to find v and n. It
%solves the sign ambiguity to have the majority of the planes in front of
%the cameras, i.e., n(:,3)>0.

function [w,v,n,output]=homographyContinuousEstimateNuclearNorm(H,E,varargin)
methodConstraint='relaxed';
lambda=0.01;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'methodconstraint'
            ivarargin=ivarargin+1;
            methodConstraint=lower(varargin{ivarargin});
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NFrames=max(E(:,1));
NPlanes=max(E(:,2));


NPairs=size(H,3);

%estimate MOuter=v(:)*n(:)'
switch lower(methodConstraint)
    case 'exact'
        idxFrames=reshape(1:3*NFrames,3,NFrames);
        idxPlanes=reshape(1:3*NPlanes,3,NPlanes);
        cvx_begin quiet
            variable MOuter(NFrames*3,NPlanes*3)
            minimize norm_nuc(MOuter)
            subject to
                for iPair=1:NPairs
                    iFrame=E(iPair,1);
                    iPlane=E(iPair,2);
                    H(:,:,iPair)+H(:,:,iPair)'==MOuter(idxFrames(:,iFrame),idxPlanes(:,iPlane))...
                        +MOuter(idxFrames(:,iFrame),idxPlanes(:,iPlane))';
                end
        cvx_end
    case 'relaxed'
        JFrames=reshape(eye(3*NFrames),[],3,NFrames);
        JPlanes=reshape(eye(3*NPlanes),[],3,NPlanes);
        cvx_begin quiet
            variable MOuter(NFrames*3,NPlanes*3)
            f=0;
            for iPair=1:NPairs
                iFrame=E(iPair,1);
                iPlane=E(iPair,2);
                M1=JFrames(:,:,iFrame)'*MOuter*JPlanes(:,:,iPlane);
                residual=H(:,:,iPair)+H(:,:,iPair)'-(M1+M1');
                f=f+residual(:)'*residual(:);
            end
            
            minimize lambda*norm_nuc(MOuter)+f; 
        cvx_end
end

%factorize MOuter into v and n
[U,S,V]=svd(MOuter);
v=S(1,1)*reshape(U(:,1),3,[]);
n=reshape(V(:,1),3,[]);

%compute sign for making sure that planes are in front of the camera
s=median(sign(n(3,:)));
if s==0
    s=1;
end

%normalize solution so that the first plane has unit distance and correct sign
kNormalization=norm(n(:,1))*s;
v=v*kNormalization;
n=n/kNormalization;

%estimate rotational velocities
w=zeros(3,NFrames);
for iFrame=1:NFrames
    flagFrame=E(:,1)==iFrame;
    idxPlane=E(flagFrame,2);
    wHat=mean(H(:,:,flagFrame)-reshape(v(:,iFrame)*vec(n(:,idxPlane))',3,3,[]),3);
    w(:,iFrame)=-vee3(wHat);
end

if nargout>2
    output.MOuter=MOuter;
    output.kNormalization=kNormalization;
    output.s=s;
end


end
