function FactorGraph_plot(x,E)
D=size(x,1);
if ~isempty(E)
    Npins=max(E(:,1));
    Nedges=size(E,1);
    hold on;
    if D==2
        scatter(x(1,1:Npins),x(2,1:Npins),150,'black','o','fill');
        scatter(x(1,Npins+1:end),x(2,Npins+1:end),200,'r','s','fill');
        grid on;
        for iEdge=1:Nedges    
            plotPoints(x(:,E(iEdge,:)),'-','Color','black')
        end
        elseif D==3
        scatter3(x(1,1:Npins),x(2,1:Npins),x(3,1:Npins),150,'black','o','fill');
        scatter3(x(1,Npins+1:end),x(2,Npins+1:end),x(3,Npins+1:end),200,'r','s','fill');
        grid on;
        for iEdge=1:Nedges    
            plotPoints(x(:,E(iEdge,:)),'-','Color','black')
        end
    end
else
    if D==2
    scatter(x(1),x(2),200,'r','s','fill');
    elseif D==3
        scatter(x(1),x(2),x(3),200,'r','s','fill');
    end
end
end