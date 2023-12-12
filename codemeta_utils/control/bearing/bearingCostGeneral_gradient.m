%Compute the gradient and Hessian of the total bearing cost
%function [gradPhi,DgradPhi]=bearingCostGeneral_gradient(y,yg,funs,ny)
%Note: if only the gradient is needed, ny can be omitted
%Inputs
%   y           [d Ny] array of current robot's bearing vectors
%   yg          [d Nyg] array of bearing vectors at the home location
%   funs        reshaping function and derivatives
%Optional arguments
%   'scales',ny range information corresponding to y (needed only for DgradPhi)
%   'matrix'    see description of outputs
%   'nosum'     do not sum contributions from all the bearings. This
%               argument is automatically implied by the 'matrix' option
%Outputs
%   gradPhi     gradient of the cost w.r.t. x. If the 'nosum' argument is
%               passed, it will have dimension [d Ny Nyg] if 
%               the optional argument 'matrix' is used (i.e., every column
%               of y is compared to every column of yg), [d Ny] if the
%               (i.e., each column of y is compared only to the
%               corresponding column of yg). If the argument 'nosum' is NOT
%               passed, the output is [d 1].
%   DgradPhi    differential of the gradient of the cost, with dimensions
%               [d d Ny Nyg] or [d d Ny] (as above).
%Note: the output are squeeze()'d, so singleton dimensions (i.e., when
%Ny==1 or Nyg==1) are removed.
%In general, in order to compute DgradPhi, the range information ny is
%necessary. However, if ny is omitted and DgradPhi is still required, an
%approximation is given by using ranges which are all equal to one.

function [gradPhi,DgradPhi]=bearingCostGeneral_gradient(y,yg,funs,varargin)
flagComputeDGrad=false;
flagMatrixOutput=false;
flagSumOutput=true;

f=funs.f;
df=funs.df;
ddf=funs.ddf;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'scales'
            ivarargin=ivarargin+1;
            ny=varargin{ivarargin};
        case 'matrix'
            flagMatrixOutput=true;
            flagSumOutput=false;
        case 'nosum'
            flagSumOutput=false;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NY=size(y,2);
NYg=size(yg,2);
d=size(y,1);

if nargout>1
    flagComputeDGrad=true;
    if ~exist('ny','var') || isempty(ny)
        ny=ones(1,NY);
    end
end



if ~flagMatrixOutput
    gradPhi=zeros(d,NY);
    if flagComputeDGrad
        DgradPhi=zeros(d,d,NY);
    end
else
    gradPhi=zeros(d,NY,NYg);        
    if flagComputeDGrad
        DgradPhi=zeros(d,d,NY,NYg);
    end
end

if ~flagMatrixOutput
    for iX=1:NY
        yi=y(:,iX);
        ygi=yg(:,iX);
        if ~flagComputeDGrad
            gradi=computePair(yi,ygi,f,df,ddf,flagComputeDGrad);
        else
            nyi=ny(iX);
            [gradi,Dgradi]=computePair(yi,ygi,f,df,ddf,flagComputeDGrad,nyi);
        end
        gradPhi(:,iX)=gradi;
        if flagComputeDGrad
            DgradPhi(:,:,iX)=Dgradi;
        end
    end
else
    for iX=1:NY
        yi=y(:,iX);
        for iXg=1:NYg
            ygi=yg(:,iX);
            if ~flagComputeDGrad
                gradi=computePair(yi,ygi,f,df,ddf,flagComputeDGrad);
            else
                nyi=ny(iX);
                [gradi,Dgradi]=computePair(yi,ygi,f,df,ddf,flagComputeDGrad,nyi);
            end
            gradPhi(:,iX,iXg)=gradi;
            if flagComputeDGrad
                DgradPhi(:,:,iX,iXg)=Dgradi;
            end
        end
    end
end

if flagSumOutput
    nsz=length(size(gradPhi));
    gradPhi=sum(gradPhi,nsz);
    if flagComputeDGrad
        DgradPhi=sum(DgradPhi,nsz+1);
    end
end
gradPhi=squeeze(gradPhi);
if flagComputeDGrad
    DgradPhi=squeeze(gradPhi);
end
end

function [gradi,Dgradi]=computePair(yi,ygi,f,df,ddf,flagComputeDGrad,nyi)
        ci=bearingComputeCosine(yi,ygi);
        fci=f(ci);
        dfci=df(ci);
        if ~flagComputeDGrad
            gradi=bearingCostGeneral_terms(yi,ygi,fci,dfci,ci);
        else
            ddfci=ddf(ci);
            [gradi,Hi]=bearingCostGeneral_terms(yi,ygi,fci,dfci,ci,ddfci);
            Dgradi=Hi/nyi;
        end
end
