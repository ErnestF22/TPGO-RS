%function A=adjrandomswitch(A,p,varargin)
%Given an adjacency matrix A, remove each link with probability p
%Optional argument:
%   connected   Remove each link only if the graph stays connected
%   symmetric   Remove undirected links (A stays symmetric)
function A=adjrandomswitch(A,p,varargin)

%rand('state',3)

flagconnected=false;    %A stays connected
flagsymmetric=false;    %A stays symmetric

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'connected'
            flagconnected=true;
        case 'symmetric'
            flagsymmetric=true;
        otherwise
            disp(varargin{ivarargin})
            error('Optional argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%find non-zero entries of A
[I,J]=find(A>0);

%for the symmetric case, consider only lower half of A
if(flagsymmetric)
    valididx=I>=J;
    I=I(valididx);
    J=J(valididx);
end

N=length(I);
K=size(A,1);

%scramble order in which to consider links
rp=randperm(N);

I=I(rp);
J=J(rp);

for(ilink=1:N)
    if(rand<p)
        if(~flagconnected)
            A(I(ilink),J(ilink))=0;
            if(flagsymmetric)
                A(J(ilink),I(ilink))=0;
            end
        else
            A1=A;
            A1(I(ilink),J(ilink))=0;
            if(flagsymmetric)
                A1(J(ilink),I(ilink))=0;
            end
            
            L=diag(sum(A1,2))-A1;       %laplacian
            
            %check if it is still connected
            if(max(components(sparse(A1)))==1)% sum(A1(I(ilink),:))>0 && sum(A1(:,J(ilink)))>0 && sum(A1(J(ilink),:))>0 && sum(A1(:,I(ilink)))>0)
                A=A1;
            end
        end
    end
end


