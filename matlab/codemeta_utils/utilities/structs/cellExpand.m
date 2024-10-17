%Expand all nested cell arrays
function c=cellExpand(c)

while any(countCellArray(c)>0)
    cExp={};
    for iCell=1:length(c)
        if ~iscell(c{iCell})
            cExp=[cExp c{iCell}];
        else
            cExp=[cExp c{iCell}{:}];
        end
    end
    c=cExp;
end
   
function count=countCellArray(c)
count=cellfun(@countCell,c);

function count=countCell(x)
if ~iscell(x)
    count=0;
else
    count=length(x);
end

function xExp=expandCell(x)
if ~iscell(x)
    xExp=x;
else
    xExp=x{:};
end