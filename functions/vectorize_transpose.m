function outmat = vectorize_transpose(inmat)
%VECTORIZE_TRANSPOSE Returns a column vector output of the vectorization
%operator after transposing input
inmat_transp = inmat';
outmat = inmat_transp(:);
end

