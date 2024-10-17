%function [E,cycles]=edgeGallery(N,type)
%Generate different patterns of edges for a graph with N nodes along with a
%cycle basis
%To indicate cycles: abs(cycle)=order of the edges, sign(cycle)=relative
%direction of cycle and edge
function [E,cycles]=edgeGallery(N,type)
switch(type)
    case 'complete'
        E=generateCompleteGraphEdges(N);
        cycles=[];
    case 'ring'
        E=[];
        for(n=1:N)
            E=[E;n mod(n,N)+1];
        end
        cycles=[1:size(E,1)]';
    case 'doublering'
        Nhalf=floor(N/2)+1;
        N1half=N-Nhalf;
        E=edgeGallery(Nhalf,'ring');
        E2=edgeGallery(N1half+2,'ring');
        E=[E;E2(2:end,:)+Nhalf-2];
        cycles=[[[1:Nhalf]';zeros(N1half+1,1)] [zeros(Nhalf-2,1);1;0;[2:N1half+2]']];
    case 'doubleringrev'
        Nhalf=floor(N/2);
        N1half=N-Nhalf;
        E=edgeGallery(Nhalf,'ring');
        E2=edgeGallery(N1half,'ring');
        E=[E;E2(2:end,:)+Nhalf-2];
        cycles=[[[1:Nhalf]';zeros(N1half-1,1)] [zeros(Nhalf-2,1);1;0;[2:N1half]']];
        cycles(:,2)=reverse_cycle(cycles(:,2));
end
