function rotationProcrustes_test
NRotations=10;

R2=rot_randn([],[],NRotations);
RRef=rot_randn();

R1Left=zeros(size(R2));
R1Right=zeros(size(R2));

for iRotation=1:NRotations
    R1Left(:,:,iRotation)=RRef*R2(:,:,iRotation);
    R1Right(:,:,iRotation)=R2(:,:,iRotation)*RRef;
end

disp([RRef;...
    rotationProcrustes(R1Left,R2,'left'); ...
    rotationProcrustes(R1Right,R2,'right'); ...
    ]);

disp(R1Left-rotationProcrustesAlign(R1Left,R2))
