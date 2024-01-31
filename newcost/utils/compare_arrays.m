function diff = compare_arrays(a,b)
%COMPARE_ARRAYS Returns max of absolute value, element-wise comparison
% of a(:) - b(:).
%Note that sizes are not checked.

diff = max(abs(a(:) - b(:)), [], "all");
fprintf("max(abs(a(:) - b(:)), [], all) %g\n", diff)

end %file function

