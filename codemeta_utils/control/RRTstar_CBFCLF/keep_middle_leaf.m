function T = keep_middle_leaf(T)
for i=1:size(T,2)
    t = [];
    idx = [];
    if ~isempty(T(i).position)
        for j = 1:size(T,2)
            flag = check_if_it_isnot_a_parent(T,j);
            if flag
                if ~isempty(T(j).position)
                    if (T(j).parent == i)
                        a = T(j).position;
                        b = T(i).position;
                        t = [t find_angle(a,b)];
                        idx = [idx j];
                    end
                end
            end
        end
    end
    if ~isempty(t)
        if size(t,2)>1
        ordered_t = sort(t);
        middle = ordered_t(1,floor(size(ordered_t,2)/2)+1);
        flag = ismember(t,middle);
        idx = idx(~flag);
        else
            idx = [];
        end
        %         [~,idx_max] = max(t);
        %         [~,idx_min] = min(t);
        %         min_max = idx([idx_max idx_min]);
        %         idx = idx(~ismember(idx,min_max));
        for k=1:size(idx,2)
            flag = check_if_it_isnot_a_parent(T,idx(k));
            if flag
                T(idx(k)).position = [];
                T(idx(k)).parent = [];
            end
        end
    end
end
end