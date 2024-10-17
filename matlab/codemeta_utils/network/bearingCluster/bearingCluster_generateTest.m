function [A,x]=bearingCluster_generateTest(testName,varargin)
testName=lower(testName);
d=2;
methodPosition='circle';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'd'
            ivarargin=ivarargin+1;
            d=varargin{ivarargin};
        case 'methodposition'
            ivarargin=ivarargin+1;
            methodPosition=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch testName
    case {'loop-3','loop-4','loop-5','loop-6'}
        NNodes=str2double(testName(6));
        A=adjgallery(NNodes,'kneigh',1);
    case 'trapezoid'
        A=adjgallery(4,'kneigh',1);
    case 'degenerate-square'
        A=adjgallery(5,'kneigh',1);
    case {'tree-2','tree-3'}
        A=blkdiag(adjgallery(3,'kneigh',2),adjgallery(4,'kneigh',2));
        A(4,3)=1;
        A(3,4)=1;
        if strcmp(testName,'tree-3')
            A=blkdiag(A,adjgallery(5,'kneigh',2));
            A(7,8)=1;
            A(8,7)=1;
        end            
    case 'tree-node-2'
        A=blkdiag(adjgallery(3,'kneigh',2),adjgallery(4,'kneigh',2));
        A(4,3)=1;
        A(3,4)=1;
        A(1,4)=1;
        A(4,1)=1;
    case {'complex-loop-5'}
        NLoop=5;
        NAdd=2;
        A=blkdiag(adjgallery(NLoop,'kneigh',1),zeros(NAdd));
        A(NLoop-1:end,NLoop-1:end)=ones(NAdd+2)-eye(NAdd+2);
    case {'complex-loop-5-2'}
        NLoop=5;
        NAdd=2;
        AComplete=ones(NAdd+2)-eye(NAdd+2);
        A=blkdiag(adjgallery(NLoop,'kneigh',1),zeros(2*NAdd));
        A(NLoop-1:NLoop+NAdd,NLoop-1:NLoop+NAdd)=AComplete;
        A([1 2 NLoop+NAdd+1:end],[1 2 NLoop+NAdd+1:end])=AComplete;
    case {'loop-rigid-4'}
        NLoop=4;
        LLoop=4;
        AComplete=ones(NLoop)-eye(NLoop);
        A=zeros(NLoop*LLoop-(LLoop-1));
        for k=1:(NLoop-1):(size(A,2)-1)
            A(k:k+NLoop-1,k:k+NLoop-1)=AComplete;
        end
        A=A(1:(end-1),1:(end-1));
        A(1,end-NLoop+2:end)=1;
        A(end-NLoop+2:end,1)=1;
    case {'butterfly','butterfly-bend'}
        A=zeros(7);
        A(1:4,1:4)=ones(4)-eye(4);
        A(1,4)=0;
        A(4,1)=0;
        A(4:7,4:7)=ones(4)-eye(4);
        A(4,7)=0;
        A(7,4)=0;
    case {'butterfly-full','butterfly-full-bend'}
        A=zeros(7);
        A(1:4,1:4)=ones(4)-eye(4);
        A(4:7,4:7)=ones(4)-eye(4);
    case {'loop-4+complete','loop-5+complete'}
        NLoop=str2num(testName(6));
        A=adjgallery(NLoop,'kneigh',1);
        NComplete=5;
        AComplete=ones(NComplete)-eye(NComplete);
        A(NLoop-1:NLoop-1+NComplete-1,NLoop-1:NLoop-1+NComplete-1)=AComplete;
    case {'bowtie'}
        AComplete=ones(4)-eye(4);
        A=blkdiag(AComplete,AComplete);
        A(4,5)=1;
        A(5,4)=1;
    case {'3-connected','3-connected-skew'}
        AComplete=ones(4)-eye(4);
        A=blkdiag(AComplete,AComplete);
        A(3,5)=1;
        A(5,3)=1;
        A(2,6)=1;
        A(6,2)=1;
        A(4,8)=1;
        A(8,4)=1;
    case 'butterfly-rigid'
        AComplete=ones(4)-eye(4);
        A=zeros(7);
        A(1:4,1:4)=AComplete;
        A(4:7,4:7)=AComplete;
        A(2,5)=1;
        A(5,2)=1;
    case 'triangle-bridge'
        AHalf=ones(4)-eye(4);
        A=blkdiag(AHalf,AHalf);
        A(2,5)=1;
        A(5,1)=1;
        A(3,5)=1;
        A(5,3)=1;
        A(8,2)=1;
        A(2,8)=1;
    otherwise
        error('Dataset name not recognized')
end

NNodes=size(A,1);

switch testName
    case {'butterfly','butterfly-rigid'}
        x=[ 0 0;
            0 1;
            1 0;
            1 1;
            1 2;
            2 1;
            2 2]';
    case 'butterfly-bend'
        d=1/sqrt(2);
        x=[ 0 0;
            0 1;
            1 0;
            1 1;
            1+d 1+d;
            1+d 1-d;
            1+2*d 1]';
    case 'trapezoid'
        x=[ 0 0;
            3 0;
            2 1;
            1 1]';
    case 'degenerate-square'
        x=[ 0 0;
            1 0;
            2 0;
            2 2;
            0 2]';
    case 'bowtie'
        x=[ 0 0;
            0 1;
            1 0;
            1 1;
            2 2;
            3 2;
            2 3;
            3 3]';
    case {'3-connected','3-connected-skew'}
        x1=[0 0;
            1 0;
            1 1;
            0 1]';
        x=[x1 x1+1.5];
        if strcmpi(testName,'3-connected-skew')
            x(1,5)=1.3;
        end
    case 'triangle-bridge'
        x1=[0 0;
            1 0;
            1 1;
            0 1]';
        x=[x1 [x1(1,:)+2; x1(2,:)]];
    otherwise
        switch lower(methodPosition)
            case 'rand'
                x=randn(d,NNodes);
            case 'circle'
                t=2*pi*(1:NNodes)/NNodes+1;
                x=[cos(t); sin(t)];
        end
end

if d>2 && size(x,1)==2
    x=[x; rand(1,NNodes)];
end
