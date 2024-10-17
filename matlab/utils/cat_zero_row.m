function A_out = cat_zero_row(A_in, num_added_rows)
%CAT_ZERO_ROW A_out is a deep copy of A_in input 2D matrix, with an extra 
% 0s row at the bottom
% One row is added by default i.e., if only one argument is passed
if (nargin < 2)
    A_out = [A_in; zeros(1, size(A_in,2))];
else
    if num_added_rows == 0
        A_out = A_in;
        return;
    end
    A_out = [A_in; zeros(num_added_rows, size(A_in,2))];
end

end %function

