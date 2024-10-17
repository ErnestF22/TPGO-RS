function POCDirectedNetworkStability
%resetRands()
NNodes=25;
%Start with a directed line graph (ensures weak connectivity)
E=[1:(NNodes-1); 2:NNodes]';
%Prepare list of all edges that can be added
EAvailable=nchoosek(1:NNodes,2);
%remove self loops
EAvailable(EAvailable(:,1)==EAvailable(:,2),:)=[];
%add both directions for edges
EAvailable=[EAvailable; fliplr(EAvailable)];
%remove edges that are already in E
EAvailable=setdiff(EAvailable,E,'rows');


while ~isempty(EAvailable)
    iEdge=randi(size(EAvailable,1));
    E=[E;EAvailable(iEdge,:)];
    EAvailable(iEdge,:)=[];
    
    [flagStable,evals]=isDirectedNetworkStable(E);
    
    if ~flagStable
        disp('Found unstable network')
        display(evals)
        keyboard
    end
    
    G=directedNetworkStabilityMatrix(E);
    GSym=directedNetworkStabilityMatrix(edgesSymmetrize(E));
    symm=@(A) (A+A')/2;
    evals=eig(symm(G'*GSym));
    flagGradient=any(real(evals)<0);
    
    if flagGradient
        disp('Found negative definite Lyapunov derivative')
        display(evals)
        keyboard
    end
    
    
end
disp('All networks passed all the tests')
