%function [Rt,vt]=rot_geodFun(R,v)
%Generate functions for the rotations Rt(t) and tangent vectors vt(t) that
%describe the geodesic rot_exp(R,t*v). If R is omitted or empty, generate a
%random R. If v is omitted, generate a random normal geodesic starting from
%R.
%
%See also rot_randGeodFun
function [Rt,dRt,R0,dR0,vVec,ddRt,dvVec]=rot_geodFun(R0,dR0,varargin)
s=1;
NGeodesics=size(dR0,3);

if ~exist('R0','var') || isempty(R0)
    R0=rot_randn(eye(3));
end
if ~exist('dR0','var') || isempty(dR0)
    dR0=rot_randTangentNormVector(R0);
end

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'speed'
            ivarargin=ivarargin+1;
            s=varargin{ivarargin};
            if ischar(s)
                switch s
                    case 'quadratic'
                        s=@(t) t^2/2;
                        ds=@(t) t;
                        dds=@(t) 1;
                    case 'cubic'
                        s=@(t) t^3/3;
                        ds=@(t) t^2;
                        dds=@(t) 2*t;
                    case 'rand'
                        s=rand;
                    otherwise
                        error('Speed profile not recognized')
                end
            end
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

u=ones(1,1,NGeodesics);
flagSpeedIsAFunction=isa(s,'function_handle');

%vectorize s if multiple geodesics are required
if ~flagSpeedIsAFunction && length(s)==1 && NGeodesics>1
    s=s*u;
elseif flagSpeedIsAFunction && length(s(0))==1 && NGeodesics>1
    s=@(t) s(t)*u;
    ds=@(t) ds(t)*u;
    dds=@(t) dds(t)*u;
end

%ensure that s, ds, and dds are function handles
if ~isa(s,'function_handle')
    dds=@(t) 0;
    ds=@(t) s;
    s=@(t) s*t; %Note: we overwrite s, so we cannot swap with line before
end

Rt=@(t) rot_exp(R0,multiprod(s(t),dR0));

vVecHat=multiprod(multitransp(R0),dR0);
%The following is the vectorized version of ds(t)*Rt(t)*R0'*dR0;
dRt=@(t) multiprod(ds(t),multiprod(Rt(t),vVecHat));
dR0=dRt(0);

if ~flagSpeedIsAFunction
    vVec=@(t) rot_vee(R0,dR0);
else
    vVec=@(t) ds(t)*rot_vee(R0,dR0);
end

ddRt=@(t) multiprod(dds(t),multiprod(Rt(t),vVecHat))...
    + multiprod(ds(t).^2,multiprod(Rt(t),multiprod(vVecHat,vVecHat)));

dvVec=@(t) dds(t)*rot_vee(R0,dR0);

function v=multiTangentComputationSingleR(Rt,R,v,t)
RtEval=Rt(t);
for iNR=1:size(v,3)
    v(:,:,iNR)=RtEval(:,:,iNR)*R'*v(:,:,iNR);
end

function v=multiTangentComputation(Rt,R,v,t)
RtEval=Rt(t);
for iNR=1:size(v,3)
    v(:,:,iNR)=RtEval(:,:,iNR)*R(:,:,iNR)'*v(:,:,iNR);
end

