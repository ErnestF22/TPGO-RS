function [Eprimeval] = PinGraph(E,x,pins,nodemembership)
%This function plots the PinGraph: Connects all of the pins associated a
%component to each other, and connects pins corresponding to neighboring
%rigid components
NPins = sum(pins);
Eprime = zeros(0,2);
Ecut = zeros(0,2);
%E = sort(E,2);
if NPins==1
    disp('SINGLE PIN GRAPH');
    Eprimeval = [];
else
    Pinidx = find(pins);
    for idxpins1 = 1:NPins 
        for idxpins2 = idxpins1+1:NPins
            temp = intersect(nodemembership{Pinidx(idxpins1)},nodemembership{Pinidx(idxpins2)});
            if ~isempty(temp)
                Eprime = [Eprime; idxpins1,idxpins2];               %New edge index for pins
                Ecut = [Eprime; Pinidx(idxpins1),Pinidx(idxpins2)]; %Orignial edge index for pins
            end
        end
    end
Eprimeval = Eprime;    
end

end