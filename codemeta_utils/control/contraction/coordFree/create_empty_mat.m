function emptyMat = create_empty_mat(flagSymbolic,n)
    if flagSymbolic
        emptyMat = sym(zeros(n));
    else
        emptyMat = zeros(n);
    end
end