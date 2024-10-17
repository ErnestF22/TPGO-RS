function POCDoubleQuadraticEmbedding
for order=2:10
    l=sym('l',[order 1]);
    tl=trifocal_line_lambda_embedding(l);
    ttl=trifocal_line_lambda_embedding(tl);
    ttlStrSorted=sortrows(sym2char(ttl));
    nConstraintsTl=sum(all(ttlStrSorted==[ttlStrSorted(2:end,:);repmat(' ',1,size(ttlStrSorted,2))],2));
    disp([order nConstraintsTl length(tl) floor((-1+sqrt(1+8*nConstraintsTl))/2)])
end
