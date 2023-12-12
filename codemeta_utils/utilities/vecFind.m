function idx=vecFind(A)
NDims=length(size(A));
if NDims==2 && size(A,2)==1
    idx=find(A);
else
    strVars=num2str(1:NDims,'I%d,');
    strVars=['[' strVars(1:end-1) ']'];
    eval([ strVars '=ind2sub(size(A),find(A));'])
    eval(['idx=' strVars ';'])
end
