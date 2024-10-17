function e=generateCompleteGraphEdges(N)
e=[];
for(i=1:N-1)
    for(j=i+1:N)
        e=[e;i j];
    end
end
