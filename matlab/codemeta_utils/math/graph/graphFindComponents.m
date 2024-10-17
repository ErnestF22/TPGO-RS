function [C]=findComponents(E,N)
if ~exist('N','var') || isempty(N)
    N=max(E(:));
end

C=1:N;

%join components and remove redundant edges
%untile there are no edges left
while ~isempty(E)
    c1=C(E(1,1));
    c2=C(E(1,2));
    
    C(C==c2)=c1;
    
    E(C(E(:,1))==c1 & C(E(:,2))==c1,:)=[];
end

%remap component numbers
L=unique(C);
C=mapValues(C,[L;1:length(L)]');
