%function A=adjgallery(N,type,opt_arg)
%Create various type of adjacency matrices for undirected graphs
%   full        completely connected, weights == 1
%   banded      fill opt_arg super and sub diagonals with weights == 1
%   kneigh      each node has 2*opt_arg neighbours
%   doublestoc  random double stocastic matrix, opt_arg is the seed for the rand
%   locrand     random geometric graph, opt_arg is the connectivity radius
%               (warning: if radius is too small, it might go in an
%               infinite loop, because it tries to always generate a
%               connected graph)
%
%   fullrand, bandedrand, kneighrand are the same as full, banded, kneigh
%   buth with random (but still symmetric) weights

%%AUTORIGHTS%%

function [A,varargout]=adjgallery(N,type,varargin)

if(N==1)
    A=0;
    return
end

flagmetropolis=false;

ivarargin=1;
A=zeros(N);
switch lower(type)
    case 'full'
        A=ones(N);
    case 'fullrand'
        A=rand(N);
        A=(A'+A)/2; % make it symmetric
    case 'banded'
        k=varargin{ivarargin};
        ivarargin=ivarargin+1;
        for(i=1:N)
            for(j=i+1:N)
                if(abs(i-j)<=k)
                    A(i,j)=1;
                    A(j,i)=1;
                end
            end
        end
    case 'tree'
        k=varargin{ivarargin};
        ivarargin=ivarargin+1;
        parents=1;  %FIFO of nodes to which we need to add children
        jnode=2;    %next available node to be assigned as a children
        while jnode<N && ~isempty(parents)
            inode=parents(1);
            parents=parents(2:end);
            
            children=min(jnode:jnode+k-1,N);
            A(inode,children)=1;
            A(children,inode)=1;
            
            parents=[parents children];
            
            jnode=jnode+k;
        end
        
        
    case 'bandedrand'
        k=varargin{ivarargin};
        ivarargin=ivarargin+1;
        for(i=1:N)
            for(j=i+1:N)
                if(abs(i-j)<=k)
                    A(i,j)=rand(1);
                    A(j,i)=A(i,j);
                end
            end
        end
    case 'kneigh'
        k=varargin{ivarargin};
        ivarargin=ivarargin+1;
        for(i=1:N)
            for(j=i+1:N)
                if(abs(i-j)<=k || abs(i-j+N)<=k)
                    A(i,j)=1;
                    A(j,i)=A(i,j);
                end
            end
        end
    case 'dkneigh'
        k=varargin{ivarargin};
        ivarargin=ivarargin+1;
        for(i=1:N)
            for(j=i+1:N)
                %if(abs(i-j)<=k || abs(i-j+N)<=k)
                if(abs(i-j)<=k)
                    A(i,j)=1;
                    %A(j,i)=A(i,j);
                end
                if(abs(i-j+N)<=k)
%                    A(i,j)=1;
                    A(j,i)=1;
                end
            end
        end
    case 'kneighrand'
        k=varargin{ivarargin};
        ivarargin=ivarargin+1;
        for(i=1:N)
            for(j=i+1:N)
                if(abs(i-j)<=k || abs(i-j+N)<=k)
                    A(i,j)=rand(1);
                    A(j,i)=A(i,j);
                end
            end
        end
    case 'sparse'
        N2=N*N;
        Nones=round(N2*varargin{1});
        order=randperm(N2);
        A=zeros(N);
        A(order(1:Nones))=1;
    case 'sparserand'
        N2=N*N;
        Nones=round(N2*varargin{1});
        order=randperm(N2);
        A=zeros(N);
        A(order(1:Nones))=rand(1:Nones,1);
    case 'doublestoc'
        if(length(varargin)>0)
            seed=varargin{1};
            rand('state',seed);
        end
        A=zeros(N);
        for(row=1:N)
            r=rand(1,N-row+1);
            A(row,row:N)=r/sum(r)*(1-sum(A(row,1:row-1)));
            A(:,row)=A(row,:);
        end
    case 'locrand'
        rratio=varargin{ivarargin};
        ivarargin=ivarargin+1;
        A=zeros(N);
        while(min(sum(A))==0)
            [A,loc]=generateGeomRandNetwork(N,rratio);
        end
        varargout{1}=loc;
    case 'otherwise'
        error(['Graph topology ''' type ''' unknown'])
end
if(~strcmp(type,'doublestoc'))
    A=A-diag(diag(A)); %zeros on the diagonal
end

%optional parameters
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'metropolis'
            flagmetropolis=true;

        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
        
if(flagmetropolis)
    d=sum(A,2);
    for(inode=1:N)
        for(jnode=1:N)
            if(A(inode,jnode)~=0)
               A(inode,jnode)=1/(1+max(d(inode),d(jnode)));
            end
        end
    end
end

function [A,loc]=generateGeomRandNetwork(N,rratio)
loc=rand(N,2);

r=rratio*sqrt(log(N)/N);

A=zeros(N);
for(inode=1:N)
    for(jnode=inode+1:N)
        if(norm(loc(inode,:)-loc(jnode,:))<r)
            A(inode,jnode)=1;
            A(jnode,inode)=1;
        end
    end
end    
    