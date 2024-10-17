function strings=sfm_poseCombine_testGetMethodString(methods)
NMethods=size(methods,1);
if NMethods>1
    strings=cell(NMethods,1);

    for iMethod=1:NMethods
        strings{iMethod}=cell2concat(cellExpand(methods{iMethod,1}));
    end
else
    strings=cell2concat(cellExpand(methodCurrent));
end
