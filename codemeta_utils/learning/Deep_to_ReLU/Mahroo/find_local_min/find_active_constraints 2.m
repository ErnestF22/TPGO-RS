function active_idx = find_active_constraints(basic,d,num_val)
% num_val is the number of total variables, d is the number of basic
% variables
val = 1:num_val; % number of total variables
active_idx = setdiff(val,basic)-d; % slack variables not in basic -> s=0 -> active
active_idx(active_idx<=0) = [];
end