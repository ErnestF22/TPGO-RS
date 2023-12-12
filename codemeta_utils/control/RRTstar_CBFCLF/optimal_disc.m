function r = optimal_disc(tree,obstacle_list,lim_x,lim_y,eta)
dimension=2;
zeta_d=pi;
miu_Xfree=freeSpace(lim_x,lim_y,obstacle_list);
gamma_rrt=2*((1+1/dimension)*(miu_Xfree/zeta_d))^(1/dimension);
card_V = size(tree,2);
r = min(gamma_rrt*((log(card_V )/card_V )^(1/3)),eta);
end