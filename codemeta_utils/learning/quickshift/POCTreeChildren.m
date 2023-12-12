function POCTreeChildren
vTree=[1;1;2;2;3;1;2;7;9;9;10];
descendents=cell(1,length(vTree));
flagContinue=true;
while flagContinue
    cntPre=cellfun(@length,descendents,'UniformOutput',true);
    descendents=appendDescendents(descendents,vTree);
    cntPost=cellfun(@length,descendents,'UniformOutput',true);
    flagContinue=any(cntPost-cntPre);
end
for id=1:length(descendents)
    fprintf('[%d] ',descendents{id});
    fprintf('\n')
end



function descendents=appendDescendents(descendents,vTree)
Nd=length(descendents);
for id=1:Nd
    descendents{vTree(id)}=union(descendents{vTree(id)},[descendents{id},id]);
end
