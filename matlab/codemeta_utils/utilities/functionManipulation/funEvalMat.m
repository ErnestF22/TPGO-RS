function r=funEvalMat(f,x,varargin)
%size of input and number of dimensions
szx=size(x);
Nszx=length(szx);

%dimension of the array corresponding to different datapoints
idxDimData=Nszx;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'indexdimensiondata'
            ivarargin=ivarargin+1;
            idxDimData=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%indeces of dimensions of a single datum
idxDimSingleDatum=setdiff(1:Nszx,idxDimData);

%number of data points
Ndata=szx(idxDimData);
%dimensions of a single datum
szSingleDatum=szx(idxDimSingleDatum);

flagReshapex=false;
if Nszx>2
    dimOrder=[idxDimSingleDatum idxDimData];
    xReshaped=reshape(permute(x,dimOrder),[],Ndata);
    flagReshapex=true;
end

%run function on the first datum in order to obtain output dimensions
flagReshapef=false;
if flagReshapex
    x0=f(reshape(xReshaped(:,1),szSingleDatum));
else
    x0=f(x(:,1));
end
szf=size(x0);

%if result is not a vector, we will need to resize at the end
%(hence we set flagReshapef=true)
Nszf=length(szf);
if Nszf>2 || szf(2)>1
    flagReshapef=true;
end

%run the function on each datum, and collect results in a 2-D array
%(even if result is multi-dim, we will take care of that afterwards)
r=zeros(prod(szf),Ndata);
for it=1:Ndata
    if flagReshapex
        xit=reshape(xReshaped(:,it),szSingleDatum);
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
    r=permute(r,dimInvOrder);
end