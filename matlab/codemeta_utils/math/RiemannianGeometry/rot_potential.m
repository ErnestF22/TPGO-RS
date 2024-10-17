function phi=rot_potential(t_node,flagmodpotential)
if(nargin<2)
    flagmodpotential=false;
end
phi=0;
for(inode=1:length(t_node))
    for(jnode=1:length(t_node))
        phi=phi+rot_dist(t_node(inode).Ri,t_node(jnode).Ri);
    end
    if(flagmodpotential)
        phi=phi+rot_dist(t_node(inode).Ri,t_node(jnode).Ri0);
    end        
end
phi=phi/2;
