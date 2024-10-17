function [Gt,vGt,G0,vG0,vGVec]=rot3r3_randGeodFun(G0,varargin)
optsG={};
s=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'speed'
            ivarargin=ivarargin+1;
            s=varargin{ivarargin};
        case 'randspeed'
            s=rand;
        case 'compact'
            optsG={'compact'};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~exist('G0','var') || isempty(G0)
    G0=RT2G(rot_randn(),randn(3,1),optsG{:});
end


N=size(G0,3);

[Rt,vt,~,v0,vVec]=rot_randGeodFun(G0(1:3,1:3,:),'speed',s);
[Tt,wt,~,w0]=real_randGeodFun(squeeze(G0(1:3,4,:)),'speed',s);
Gt=@(t) RT2G(Rt(t),Tt(t),optsG{:});
vGt=@(t) [vt(t) permute(wt(t),[1 3 2]); zeros(1,4,N)];
vG0=[v0 permute(w0,[1 3 2]); zeros(1,4,N)];
vGVec=[vVec;w0];
