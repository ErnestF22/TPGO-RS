%function g=gradG(t_node,Rmean)
%Compute the gradient of the sum of squared distances on SO(3) at Rmean
function g=gradG(t_node,Rmean)
g=0;
for(inode=1:length(t_node))
    g=g+logrot(t_node(inode).Ri'*Rmean);
end

