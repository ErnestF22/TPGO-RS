function result = check_tableau_after_flipping(z,bit,P,As,Ar,Ab_set,active_idx)
% Arf tableau from function
Arf = find_tableau_after_flip(z,bit,P,As,Ar,Ab_set);

% Arb tableau from dual simplex method
z(bit,:)=~z(bit,:); % z for another region
[Asf,basic,d] = getdualA1(z,Ab_set);
[basic,result,Arb,~]=dual_simplex1(Asf,basic,d);

num_val = length(result);
active_idxs = find_active_constraints(basic,d,num_val); % vertex from dual simplex(may not be the vertex we want)

if active_idxs == active_idx % if same vertex
    % check
    if sum(Arf-Arb)<=1e-9
        result = true;
    else
        result = false;
    end
else % if not same vertex
    [T,~,zs]=find_vertices1(basic,Arb,d);
    basic_idx = find(~sum(zs-active_idx,2)); % find the same vertex
    if sum(Arf-T{basic_idx})<=1e-9
        result = true;
    else
        result = false;
    end
end
