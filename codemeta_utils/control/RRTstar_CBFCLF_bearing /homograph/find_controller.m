function [K,k_added] = find_controller(exitD,Wb,Wl,Ah,Ax,bx,bh,xe,y,L,Cl,Cb,E,flag)
k_added = [0;0];
if flag == 0
    [K,~,~,~,~,k_added] = optFirstorderWithU(Wb,Wl,Cb,Cl,exitD,Ah,Ax,bx,bh,L,xe);
%     [K, k_added] = optFirstorderWithU_center(Wb,Wl,Cb,Ah,Ax,bx,bh,L,xe);
else
    flag_test = 1;
    [sMin,sMax,zMax,zMin] = MinMaxDis(y,L,flag_test);
%     S = ones(4,1);
%     K = optFirstorder_two_inner_opt(Wb,Wl,Cb,Cl,exitD,Ah,Ax,bx,bh,L,xe,zMax,E);
%     K =  optFirstorderWithS(Wb,Wl,Cb,Cl,exitD,Ah,Ax,bx,bh,L,xe,sMin,sMax);
%     K = optFirstorder_two_inner_opt_withoutPb(Wb,Wl,Cb,Cl,exitD,Ah,Ax,bx,bh,L,xe,Zmax,E);
%     K = optFirstorder_yimg(Wb,Wl,Cb,Cl,Ah,Ax,bx,bh,L,xe,E,exitD,zMax);
%     [K,k_added] =optFirstorder_yimg_predual(Wb,Wl,Cb,Cl,Ah,Ax,bx,bh,L,xe,E,exitD,zMax,zMin);
    [K,k_added] = optFirstorderWithS_test_inequality_trace(Wb,Wl,Cb,Cl,exitD,Ah,Ax,bx,bh,L,xe,sMin,sMax);    
end
end