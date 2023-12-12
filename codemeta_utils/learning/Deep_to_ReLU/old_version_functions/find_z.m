function idx = find_z(set_out,z)
for i=1:size(set_out,2)
    if set_out(i).z==z
        idx = i;
    end
end