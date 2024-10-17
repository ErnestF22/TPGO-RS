function A=trifocal_line_lambda_buildEquations(c,idx,V)
NLambda=size(V,2);
NT=NLambda*(NLambda+1)/2;
NConstraints=size(c,3);
LConstraints=size(c,2);
A=zeros(NConstraints,NT);
for iConstraint=1:NConstraints
    a=zeros(NT,1);
    for jConstraint=1:LConstraints
        v1=V(idx(1,jConstraint,iConstraint),:);
        v2=V(idx(2,jConstraint,iConstraint),:);
        a=a+c(1,jConstraint,iConstraint)*trifocal_line_lambda_symprod(v1,v2);
    end
    A(iConstraint,:)=a';
end
