function [min_distance,act_index] = find_flip_index(y_result,ACT)
% find the edge which contains the minimum and second minimum distance vertices
flagfind = false;
[min_distance, y_index1] = min(y_result); % find minimum distance
act1 = ACT(y_index1,:); % active constraints at this vertex
y_result(y_index1) = Inf; % checked, set it to Inf
while ~flagfind
    [~, y_index2] = min(y_result);
    act2 = ACT(y_index2,:); % active constraints at this vertex
    y_result(y_index2) = Inf; % checked, set it to Inf
    act_index = intersect(act1,act2); % intersection of active constraints of two vertices
    if ~isempty(act_index)
        flagfind = true; % find z to flip if there is intersection
    end
    % if there is no intersection, find the third minimum distance and so
    % on
    % another way is to try to flip d active constraints, we can get d
    % regions, find the closest point among d regions
end

end