function R = make_random_sym_blockrot_matrix(d, N)
%Make random block symmetric matrix, with (i,j) block equal to the
%inverse of (j,i) block; rotations are in SO(d)
    
    upper_triang_size_nodiag = N*(N-1)/2;
    R_nosym = randrot_som(d, upper_triang_size_nodiag); %NOTE: NOT all of this will be transformed into rotation matrix
    R = eye(d*N);
    i=1;
    j=2;
    max_to_reach = N;
    % fprintf("size (n-1)n/2 %g\n", size(R_nosym, 1));
    for totid = 1:size(R_nosym, 3)
        %     R((i-1)*d+1:(i  -1)*d+d, (j-1)*d+1:(j-1)*d+d) = quat2rotm(R_nosym(totid));
        %     R((j-1)*d+1:(j-1)*d+d, (i-1)*d+1:(i-1)*d+d) = inv(quat2rotm(R_nosym(totid))); % ' instead of inv() would do the same
        R((i-1)*d+1:(i-1)*d+d, (j-1)*d+1:(j-1)*d+d) = R_nosym(:,:,totid);
        R((j-1)*d+1:(j-1)*d+d, (i-1)*d+1:(i-1)*d+d) = inv(R_nosym(:,:,totid)); % ' instead of inv() would do the same
        j = j+1;
        if (j>max_to_reach)
            i = i+1;
            j = i+1;
        end
    end
    % disp(R)
end