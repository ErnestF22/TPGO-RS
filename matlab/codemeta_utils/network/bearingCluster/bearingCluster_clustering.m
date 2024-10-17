function [membership,L,s]=bearingCluster_clustering(E,u,varargin)

global cNext cInitial
tol=1e-12; 
flagSeparateComponents=true;
flagFirstCall=false;
Threshold=0.3;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        case 'flagfirstcall'
            ivarargin=ivarargin+1;
            flagFirstCall=varargin{ivarargin};
            idxvarargin=ivarargin;
        case 'flagseparatecomponents'
            ivarargin=ivarargin+1;
            flagSeparateComponents=varargin{ivarargin};
        case 'threshold'
            ivarargin=ivarargin+1;
            Threshold=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%This part clusters edges into interdependent sets
CUnsigned=grCycleBasisMultiComp(E);
C=grOrientCycleBasis(CUnsigned,E)';
M=bearingCluster_measurementMatrixFromC(C,u);
[L,s]=bearingCluster_nullSpaceBasis(M,tol);
N=bearingCluster_clusterVectors(L);
%[membership,info]=projective_quickshift(N,'threshold',0.5);
[membership,info]=projective_quickshift(N,'optsbreaktree',{'Threshold',0.3});

if flagFirstCall
    cNext=max(membership)+1;
    varargin{idxvarargin}=false;
    cInitial=1;
end

%Check if a set of interdependent edges are still interdependent with
%other edges removed
c=unique(membership);
if length(c)>1 && flagSeparateComponents
    if ~flagFirstCall
        membership=relabelSubComponents1(membership,cInitial);
    end
    c=unique(membership);
    Nc=length(c);
    for ic=1:Nc
        idxEc=find(membership==c(ic));
        Ec=E(idxEc,:);
        uc=u(:,idxEc);
        cInitial=c(ic);
        membershipEc=bearingCluster_clustering(Ec,uc,varargin{:});
        membership=relabelSubComponents2(membership,idxEc,membershipEc);
    end
else
    membership=cInitial*membership;
end
end
%--------------------------------------------------------------------------
function membershipVal=relabelSubComponents1(membership,cInitial)
global cNext
newComps=find(membership~=1);
membership(newComps)=membership(newComps)-2+cNext;
membership(membership==1)=cInitial;
cNext=max(membership)+1;
membershipVal=membership;
end
%--------------------------------------------------------------------------
function membershipVal=relabelSubComponents2(membership,idxEc,membershipEc)
%global cNext
membership(idxEc)=membershipEc;
membershipVal=membership;
% c=unique(membershipEc);
% if length(c)==1
%     membershipVal=membership;
% else
%     initMemb=unique(membership(idxEc));
%     oldComps=find(membershipEc==1);
%     newComps=find(membershipEc~=1);
%     membershipEc(newComps)=membershipEc(newComps)-2+cNext;
%     membershipEc(oldComps)=ones(1,length(oldComps))*initMemb;
%     membership(idxEc)=membershipEc;
%     cNext=max(membership)+1;
%     membershipVal=membership;
% end
end


%%%% OLD version from Roberto
%function [membership,L,s]=bearingCluster_clustering(E,u,varargin)
%tol=1e-12;
%flagSeparateComponents=true;
%
%%optional parameters
%ivarargin=1;
%while(ivarargin<=length(varargin))
%    switch(lower(varargin{ivarargin}))
%        case 'tol'
%            ivarargin=ivarargin+1;
%            tol=varargin{ivarargin};
%        case 'flagseparatecomponents'
%            ivarargin=ivarargin+1;
%            flagSeparateComponents=varargin{ivarargin};
%        otherwise
%            error(['Argument ' varargin{ivarargin} ' not valid!'])
%    end
%    ivarargin=ivarargin+1;
%end
%
%CUnsigned=grCycleBasisMultiComp(E);
%C=grOrientCycleBasis(CUnsigned,E)';
%M=bearingCluster_measurementMatrixFromC(C,u);
%
%[L,s]=bearingCluster_nullSpaceBasis(M,tol);
%%disp(s)
%
%N=bearingCluster_clusterVectors(L);
%
%[membership,info]=projective_quickshift(N,'threshold',0.5);
%
%if flagSeparateComponents
%    c=unique(membership);
%    cNext=max(c)+1;
%    Nc=length(c);
%    for ic=1:Nc
%        idxEc=find(membership==c(ic));
%        Ec=E(idxEc,:);
%        membershipEc=grCompEdges(Ec);
%        NcEc=max(membershipEc);
%        for icEc=2:NcEc
%            %relabel the different connected components in the cluster
%            membership(idxEc(membershipEc==icEc))=cNext;
%            cNext=cNext+1;
%        end
%    end
%end
