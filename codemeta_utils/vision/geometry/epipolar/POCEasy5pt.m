function POCEasy5pt
% [G1,G2,x1,x2]=epipolar_dataset();
% NPoints=5;
% x1=homogeneous(x1(:,1:NPoints),3);
% x2=homogeneous(x2(:,1:NPoints),3);
% 
% A=[];
% for iPoint=1:NPoints
%     A=[A;kron(x2(:,iPoint),x1(:,iPoint))'];
% end
% ETruth=epipolarBuildEFromG(G1,G2,'references');
% 
% %disp(A*ETruth(:))
% 
% [U,S,V]=svd(A);

% E0=V(:,end);
% E1=V(:,end-1);
% E2=V(:,end-2);
% E3=V(:,end-3);

bE=sym('bE',[9,4]);
E0=bE(:,1);
E1=bE(:,2);
E2=bE(:,3);
E3=bE(:,4);
xy=sym('xy',[2,1]);
syms z
E=reshape(E0*xy(1)+E1*xy(2)+E2*z+E3,3,3);

p=[
    collect(reshape(2*(E*E.')*E-trace(E*E.')*E,[],1),xy);
    det(E)
    ];

%disp(p)
%matlabFunction(p,'file',[mfilename '_p']);
C=[];
T=[];
for ip=1:length(p)
    [C1,T1]=coeffs(p(ip),xy);
    C=[C;C1];
    T=[T;T1];
end

% Cz=cell(10);
% Tz=cell(10);
% for iz=1:10
%     for jz=1:10
%         [Cz{iz,jz},Tz{iz,jz}]=coeffs(C(iz,jz),z);
%     end
% end
%pz=det(C);

Cz=[];
CzE=[];
for kz=1:4
    Cz=cat(3,Cz,sym(['Cz' num2str(kz) '_'],[10 10]));
    CzE=cat(3,CzE,sym(['CzE' num2str(kz) '_'],[10 10]));
end

for iz=1:10
    for jz=1:10
        c=coeffs(C(iz,jz),z);
        switch jz
            case {1,2,4,7}
                c=[0 0 0 c];
            case {3,5,8}
                c=[0 0 c];
            case {6,9}
                c=[0 c];
            case 10
                %nothing to do
            otherwise
                error;
        end
        for kz=1:4
            if isequal(c(kz),0)
                Cz(iz,jz,kz)=0;
                CzE(iz,jz,kz)=0;
            else
                CzE(iz,jz,kz)=c(kz);
            end
        end
    end
end

CCz=Cz(:,:,1)*z^3+Cz(:,:,2)*z^2+Cz(:,:,3)*z+Cz(:,:,4);

pz=detLeibniz(CCz);

save([mfilename '_data'])

