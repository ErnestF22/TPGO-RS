function r=funEvalVec(f,x,varargin)
%size of input and number of dimensions
szx=size(x);
Nszx=length(szx);

%dimension of the array corresponding to different datapoints
idxDimVector=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'indexdimensionvector'
            ivarargin=ivarargin+1;
            idxDimVector=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%index of dimension of data
idxDimData=setdiff(1:Nszx,idxDimVector);

%dimensions of data
szDimData=szx(idxDimData);

%number of data points
Ndata=prod(szDimData);

flagReshapex=false;
if Nszx>2 || idxDimVector~=1
    dimOrder=[idxDimVector idxDimData];
    xReshaped=reshape(permute(x,dimOrder),[],Ndata);
    flagReshapex=true;
end

%run function on the first datum in order to obtain output dimensions
flagReshapef=false;
if flagReshapex
    x0=f(xReshaped(:,1));
else
    x0=f(x(:,1));
end
szf=size(x0);
if length(szf)==2 && szf(2)==1
    szf=szf(1);
end

%if result is not a vector, we will need to resize at the end
%(hence we set flagReshapef=true)
if length(szf)>1
    flagReshapef=true;
end

%run the function on each datum, and collect results in a 2-D array
%(even if result is multi-dim, we will take care of that afterwards)
r=zeros(prod(szf),Ndata);
for it=1:Ndata
    if flagReshapex
        xit=xReshaped(:,it);
    else
        xit=x(:,it);
    end
    r(:,it)=reshape(f(xit),[],1);
end

%reshape result if necessary
if flagReshapef
    r=reshape(r,[szf Ndata]);        
end

%permute result if necessary to match the original ordering of x
if flagReshapex
    dimInvOrder(dimOrder)=1:length(dimOrder);
    r=permute(reshape(r,[szf szx(dimOrder(2:end))]),dimInvOrder);
end
