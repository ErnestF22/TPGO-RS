function getSeqPerm_test
method='sjc';
N=7;
h=getSeqPerm(1:N,'method',method);
NN=factorial(N);
S=zeros(NN,N);
for iS=1:NN
    s=h();
    if iS>1 && any(all(S(1:iS-1,:)==ones(iS-1,1)*s,2))
        error('Test failed')
    else
        S(iS,:)=s;
    end
end
disp('Test succeeded')

