function c1c2=kroncat(c1,c2)
if ~iscell(c1)
    c1={c1};
end
if ~iscell(c2)
    c2={c2};
end
nb_c1=numel(c1);
nb_c2=numel(c2);
c1c2=cell(nb_c1,nb_c2);

for iCell1=1:nb_c1
    for iCell2=1:nb_c2
        c1c2{iCell1,iCell2}=[c1{iCell1} c2{iCell2}];
    end
end

