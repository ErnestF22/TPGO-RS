function P=RTK2P(R,T,K,varargin)
methodAbsolutePoses='pose';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
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

if strcmpi(methodAbsolutePoses,'reference')
    [R,T]=invRT(R,T);
end
    
[d,NP]=size(T);

if ~exist('K','var')
    K=eye(d);
end

if size(K,3)<NP
    K=repmat(K,[1 1 NP]);
end

P=zeros(d,d+1,NP);
for iP=1:NP
    P(:,:,iP)=K(:,:,iP)*[R(:,:,iP) T(:,iP)];
end
