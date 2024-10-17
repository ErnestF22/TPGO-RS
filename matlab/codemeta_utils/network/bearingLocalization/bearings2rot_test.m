function bearings2rot_test
%We use reference interpretation, i.e., body-to-absolute
%ijk translations
T=[zeros(3,1),randn(3,2)];
Ti=T(:,1);
Tj=T(:,2);
Tk=T(:,3);

%ij rotations
R=rot_randn(eye(3),[],2);
Ri=eye(3);%R(:,:,1);
Rj=R(:,:,2);

%relative rotation/translations
[Rij,Tij]=computeRelativePoseFromRT(Ri,Ti,Rj,Tj,'references');
%Tij is Tj expressed in Bi
[~,Tji]=computeRelativePoseFromRT(Rj,Tj,Ri,Ti,'references');
Tik=rigidTransform(Ri,Ti,Tk,'references');
Tjk=rigidTransform(Rj,Tj,Tk,'references');

%bearings
tij=cnormalize(Tij);
tji=cnormalize(Tji);
tik=cnormalize(Tik);
tjk=cnormalize(Tjk);

RijEst=bearings2rot(tij,tji,tik,tjk);

disp('Error between estimated and real relative rotations should be zero')
disp(norm(Rij-RijEst,'fro'))
