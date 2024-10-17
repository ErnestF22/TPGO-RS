function t=trifocal_line_lambda_embedding(lambda)
NLambda=length(lambda);
NT=NLambda*(NLambda+1)/2;
cnt=1;

% %check that lambda is not symbolic
% if isempty(symvar(lambda))
%    t=zeros(NT,1);
% end
for iLambda=1:NLambda
    for jLambda=iLambda:NLambda
        t(cnt)=lambda(iLambda)*lambda(jLambda);
        cnt=cnt+1;
    end
end
