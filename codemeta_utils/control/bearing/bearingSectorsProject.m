function gProj=bearingSectorsProject(g,y,a)
uProj=bearingSectorsFind(g,y,a);
if isempty(uProj) || g'*uProj<0
    gProj=zeros(2,1);
else
    gProj=g'*uProj*uProj;
end
