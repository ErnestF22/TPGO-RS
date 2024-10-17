function vw=trifocal_line_lambda_symprod(v,w)
NLambda=length(v);
NT=NLambda*(NLambda+1)/2;
cnt=1;
vw=zeros(NT,1);
for iLambda=1:NLambda
    for jLambda=iLambda:NLambda
        if iLambda==jLambda
            vw(cnt)=v(iLambda)*w(jLambda);
        else
            vw(cnt)=v(iLambda)*w(jLambda)+v(jLambda)*w(iLambda);
        end
        cnt=cnt+1;
    end
end
