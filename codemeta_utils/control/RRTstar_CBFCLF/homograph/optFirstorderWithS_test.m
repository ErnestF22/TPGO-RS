%this function finds K1
function K = optFirstorderWithS_test(Wb,Wl,Cb,Cl,z,Ah,Ax,bx,bh,L,xe,S,sMin,sMax)
%Wb:steering factor for cbf
%Wl:steering factor for clf
%Cb:cbf coefficoient
%Cl:clf coefficoient
%z:exit direction
%Ah:barrier matrix
%bh:barrier bias
%Au:velocity matrix
%bu:velocity bias
%Ax:convex cell matrix
%bx:convex cell bias
%Y: vectorized of landmarks
%d:dimestion

I = repmat(eye(2),4,1);
nBarriers = size(Ah,1);
nConvex = size(Ax,1);
L1 = blkdiag(L(:,1),L(:,2),L(:,3),L(:,4));
L2 = reshape(L,[],1);
nLandmark = length(L);

cvx_begin

variables  delta_clf delta_b K(2,nLandmark*2)...
           delta_eqCBF(2,3) delta_eqCLF(2)...
           lambda_cbf(nConvex,nBarriers) lambda_clf(nConvex,1) ...
           P1(4,3) P2(4,3) P3(4,3) P4(4,3) P5(4,3) P6(4,3)... 
           P7(4,3) P8(4,3) P9(4,3) P10(4,3) ...
           P11(4,1) P12(4,1) P13(4,1) P14(4,1) P15(4,1) P16(4,1)...
           P17(4,1) P18(4,1) P19(4,1) P20(4,1)
       
Kbx = [K(:,1) K(:,3) K(:,5) K(:,7)];
Kby = [K(:,2) K(:,4) K(:,6) K(:,8)];

delta_cbf = [delta_b;delta_b;delta_b];

minimize(delta_b+delta_clf+sum(sum(delta_eqCBF))+sum(delta_eqCLF))

subject to
%%
%constraints for safety:
% max_S: 
% bx'*lambda_cbf <= (delta_cbf+Ah*K*S*L2+Cb*bh)'
% s.t: smin<=S<=smax
[-sMin' sMax']*[P1;P2] <= (delta_cbf+Cb*bh)'-bx'*lambda_cbf
[-eye(nLandmark) eye(nLandmark)]*[P1;P2]>=(-Ah*K*L1)'
P1>=0
P2>=0
%%
% Ax'*lambda_cbf == (Ah*K*S*I-Ah*Cb)'
% % % split into two inequalities:

%  Ah*Kbx*ones(4,1) <=  ([1 0]*((Ah*Cb)'+Ax'*lambda_cbf))'+delta_eqCBF(1,:)'
% -Ah*Kbx*ones(4,1) <= -([1 0]*((Ah*Cb)'+Ax'*lambda_cbf))'+delta_eqCBF(1,:)'
% 
%  Ah*Kby*ones(4,1) <=  ([0 1]*((Ah*Cb)'+Ax'*lambda_cbf))'+delta_eqCBF(2,:)'
% -Ah*Kby*ones(4,1) <= -([0 1]*((Ah*Cb)'+Ax'*lambda_cbf))'+delta_eqCBF(2,:)'
 %%
[-sMin' sMax']*[P3;P4] <= [1 0]*((Ah*Cb)'+Ax'*lambda_cbf)+delta_eqCBF(1,:)
[-eye(nLandmark) eye(nLandmark)]*[P3;P4]>= (Ah*Kbx)'

[-sMin' sMax']*[P5;P6] <= -[1 0]*((Ah*Cb)'+Ax'*lambda_cbf)+delta_eqCBF(1,:)
[-eye(nLandmark) eye(nLandmark)]*[P5;P6]>= -(Ah*Kbx)'

P3>=0
P4>=0
P5>=0
P6>=0
% % 
[-sMin' sMax']*[P7;P8] <= [0 1]*((Ah*Cb)'+Ax'*lambda_cbf)+delta_eqCBF(2,:)
[-eye(nLandmark) eye(nLandmark)]*[P7;P8]>= (Ah*Kby)'

[-sMin' sMax']*[P9;P10] <= -[0 1]*((Ah*Cb)'+Ax'*lambda_cbf)+delta_eqCBF(2,:)
[-eye(nLandmark) eye(nLandmark)]*[P9;P10]>= -(Ah*Kby)'

P7>=0
P8>=0
P9>=0
P10>=0

%%
% % constraints for stability:
% maax_S:
% bx'*lambda_clf<= (delta_clf-z'*K*S*L2+Cl*z'*xe)'
% s.t: smin<=S<=smax
[-sMin' sMax']*[P11;P12]<= (delta_clf+Cl*z'*xe)'-bx'*lambda_clf
[-eye(nLandmark) eye(nLandmark)]*[P11;P12]>=(z'*K*L1)'
P11>=0
P12>=0

%%
% % % Ax'*lambda_clf == (-z'*K*S*I+z'*Cl)'
% split into two inequalities:
%  z'*Kbx*ones(4,1) <=  [1 0]*((z'*Cl)'-(Ax'*lambda_clf))+delta_eqCLF(1)
% -z'*Kbx*ones(4,1) <= -[1 0]*((z'*Cl)'-(Ax'*lambda_clf))-delta_eqCLF(1)
% 
%  z'*Kby*ones(4,1) <=  [0 1]*((z'*Cl)'-(Ax'*lambda_clf))+delta_eqCLF(2)
% -z'*Kby*ones(4,1) <= -[0 1]*((z'*Cl)'-(Ax'*lambda_clf))-delta_eqCLF(2)
%%
[-sMin' sMax']*[P13;P14] <= [1 0]*((z'*Cl)'-(Ax'*lambda_clf))+delta_eqCLF(1)
[-eye(nLandmark) eye(nLandmark)]*[P13;P14]>=(z'*Kbx)'

[-sMin' sMax']*[P15;P16] <= -[1 0]*((z'*Cl)'-(Ax'*lambda_clf))+delta_eqCLF(1)
[-eye(nLandmark) eye(nLandmark)]*[P15;P16]>=-(z'*Kbx)'

P13>=0
P14>=0
P15>=0
P16>=0

[-sMin' sMax']*[P17;P18] <= [0 1]*((z'*Cl)'-(Ax'*lambda_clf))+delta_eqCLF(2)
[-eye(nLandmark) eye(nLandmark)]*[P17;P18]>=(z'*Kby)'

[-sMin' sMax']*[P19;P20] <= -[0 1]*((z'*Cl)'-(Ax'*lambda_clf))+delta_eqCLF(2)
[-eye(nLandmark) eye(nLandmark)]*[P19;P20]>=-(z'*Kby)'

P17>=0
P18>=0
P19>=0
P20>=0

%%
lambda_cbf >= 0
lambda_clf >= 0

delta_b <= 0
delta_clf <= 0
delta_b >= -0.1
delta_clf >= -0.1
%
% delta_eqCBF <= 0.0001
% delta_eqCLF <= 0.0001


cvx_end
end

