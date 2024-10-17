function A_out = cat_zero_col(A_in, num_added_cols)
%CAT_ZERO_ROW A_out is a deep copy of A_in input 2D matrix, with an extra 
% 0s row at the bottom
% One col is added by default i.e., if only one argument is passed
if (nargin < 2)
    A_out = [A_in, zeros(size(A_in,1), 1)];
else
    A_out = [A_in, zeros(size(A_in,1), num_added_cols)];
end
end %function

