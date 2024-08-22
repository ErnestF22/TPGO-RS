function [is_equal_bool, err] = is_equal_floats(A,B, thr)

if ~exist('thr', 'var')
    thr = 1e-5;
end

err = max(abs(A(:) - B(:)), [], "all");
is_equal_bool = err < thr;


end %file function
