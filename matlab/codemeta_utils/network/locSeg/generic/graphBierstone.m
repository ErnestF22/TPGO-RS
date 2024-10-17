function C=bierstone(A)
NNodes=size(A,1);

%initialization phase
flagInitialNodeFound=false;
jStart=NNodes;
while ~flagInitialNodeFound
    jStart=jStart-1;
    W=Mj(A,jStart);
    if ~isempty(W)
        flagInitialNodeFound=true;
    end
end

%init connected components
C=initComponent(jStart,W);

%grow and merge connected components
for j=jStart-1:-1:1
    M=Mj(A,j);
    if ~isempty(M)
        W=M;
        L=length(C);
        %start search in existing subgraphs
        for n=1:length(C)
            T=intersect(M,C{n});
            if length(T)>1
                W=setdiff(W,T);
                if setequal(T,C{n})
                    %C{n} and T are the same
                    %step 15
                    C{n}=[j C{n}];
                    C=cleanC(C,C{n},n+1);
                    %step 16
                    if setequal(T,M)
                        break
                    end
                else
                    %step 12
                    S=[j T];
                    if setequal(T,M)
                        C=cleanC(C,S,L+1);
                        C{end+1}=S;
                    else
                        %step 13
                        if ~isSubset(C,S)
                            %set 14
                            C=findAndReplaceOrAdd(C,S);
                        end
                    end
                end
            end
        end
        if ~isempty(W)
            %add whatever remains of W as new components
            Cnew=initComponent(j,W);
            C=[C Cnew];
        end
    end
end

%returns true if sets a and b are equal
function flag=setequal(a,b)
flag=(length(a)==length(b)) && all(sort(a)==sort(b));

%finds the first set in C which is a subset of S, replace it, and remove
%all subsequent subsets. If not found, just add at the end
function C=findAndReplaceOrAdd(C,S)
flagFound=false;
for iC=1:length(C)
    if all(ismember(C{iC},S))
        flagFound=true;
        C{iC}=S;
        C=cleanC(C,S,iC+1);
    end
end
if ~flagFound
    C{end+1}=S;
end

%returns true if T is a subset of any set in C
function flag=isSubset(C,S)
flag=false;
for iC=1:length(C)
    if all(ismember(S,C{iC}))
        flag=true;
        break
    end
end

%init components with current nodes and nodes in a given set
function C=initComponent(j,W)
NW=length(W);
C=cell(1,NW);
for iW=1:NW
    k=W(iW);
    C{iW}=[j k];
end

%set of nodes k>j connected to node j
function W=Mj(A,j)    
W=find(A(j,j+1:end))+j;

%removes all elements C{n},...,C{end} which are subsets of S
function C=cleanC(C,S,n)
N=length(C);
flagKeep=true(1,N);
for iC=n:N
    if all(ismember(C{iC},S))
        flagKeep(iC)=false;
    end
end
C=C(flagKeep);
