function a=cell2concat(c)
c=cellfun(@stringify,c,'UniformOutput',false);
a=[c{:}];

function s=stringify(a)
if ischar(a)
    s=a;
else
    if isnumeric(a) || islogical(a)
        s=num2str(a,'%d');
    end
end
