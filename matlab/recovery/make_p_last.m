function P_last_out = make_p_last(p,d)
%MAKE_P_LAST Return matrix that extracts the last p-d rows of matrix A
% when multiplied on left side of A i.e., P_last_out * A
P_last_out = [zeros(p-d,d), eye(p-d)];
end

