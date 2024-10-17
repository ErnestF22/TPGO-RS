%function [Gtruth,X]=testNetworkCreateAbsolutePoses(N, varargin)
%Generates absolute poses for N cameras (default N=7). If output argument X is requested,
%generates also 3D points
%Optional arguments
%   'NPoints'       number of 3D points to generate (default NPoints=30)
%

%%AUTORIGHTS%%

function [Gtruth,X]=testNetworkCreateAbsolutePoses(N, varargin)
if ~exist('N','var')
    N=7;
end
NPoints=30;
flagIdentityRot=false;
flagTranslPerturbation=false;
methodAbsolutePoses='pose';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'npoints'
            ivarargin=ivarargin+1;
            NPoints=varargin{ivarargin};
        case 'identityrot'
            flagIdentityRot=true;
        case 'translperturbation'
            flagTranslPertubation=true;
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

Gtruth=zeros(4,4,N);
for iCamera=1:N
    Gtruth(1:3,1:3,iCamera)=abspose2rot(0,pi/8*mod(iCamera-1,2)+pi/2,(iCamera-1)/N*2*pi);
    T=[0;0;-8];
    if flagTranslPerturbation
        T=T+2*randn(3,1);
    end
    Gtruth(1:3,4,iCamera)=Gtruth(1:3,1:3,iCamera)*T;
    if flagIdentityRot
        Gtruth(1:3,1:3,iCamera)=eye(3);
    end
    Gtruth(4,4,iCamera)=1;
    Gtruth(:,:,iCamera)=invg(Gtruth(:,:,iCamera));
end

if nargout>1
    X=9*(rand(3,NPoints)-0.5);
end

switch methodAbsolutePoses
    case 'pose'
        %nothing to do
    case 'reference'
        Gtruth=invg(Gtruth);
    otherwise
        error('methodAbsolutePoses not recognized')
end

% impixel=1000;
% %     for iCamera=1:N
% %         ximage=persp_project(dehom(t_node(iCamera).gitruth*[X; ones(1,size(X,2))]));
% %         t_node(iCamera).ximage=ximage+sigmanoisePoints/impixel*[randn(2,size(ximage,2));zeros(1,size(ximage,2))];
% %     end
% 
% ximage=zeros(3,NPoints,N);
% XTransf=zeros(3,NPoints,N);
% for iCamera=1:N
%     XTransf(:,:,iCamera)=dehom(Gitruth(:,:,iCamera)*[X; ones(1,size(X,2))]);
%     ximage(:,:,iCamera)=persp_project(XTransf(:,:,iCamera));
% end
% 
% figure(1)
% for iCamera=1:N
%     draw3dcameraFromPose(Gitruth(1:3,1:3,iCamera),Gitruth(1:3,4,iCamera));
%     hold on
% end
% plot3(X(1,:),X(2,:),X(3,:),'k.')
% hold off
% 
% figure(2)
% for iCamera=1:N
%     subplot(3,ceil(N/3),iCamera)
%     plot(ximage(1,:,iCamera),ximage(2,:,iCamera),'.');
% end
% 
% for iCamera=1:2
%     for jCamera=1:1
%         if (iCamera~=jCamera)
%             x1=ximage(:,:,iCamera);
%             x2=ximage(:,:,jCamera);
%             [R,p] = discrete_motion(x1,x2);
%             gij=Gitruth(:,:,iCamera)*invg(Gitruth(:,:,jCamera));
%             disp(rot_dist(R,gij(1:3,1:3))*180/pi)
%         end
%     end
% end
 