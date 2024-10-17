function POCTwistAmbiguity
Q=zeros(6,3);
Q(1:3,:)=rot(3*pi/8*[0;1;0]);
Q(4:6,:)=rot(5*pi/8*[0;1;0]);
X=[2;0;1.5];

for k=1:4
    QFlip=essential_flipAmbiguity(Q,k);
    
    [R1,T1,R2,T2]=essential_getRTPair(QFlip);
    
    subplot(2,4,2*k-1)    
    showPair(R1,T1,R2,T2,[num2str(k) 'p'])
    subplot(2,4,2*k)
    showPair(R1,T1,R2,-T2,[num2str(k) 'n'])
end
note=char(...
    'Equivalences:',...
    '1p==4n',...
    '1n==4p',...
    '2p==3n',...
    '2n==3p');
disp(note)


function showPair(R1,T1,R2,T2,titleText)
optsDrawCamera={'references','scale',0.3};
optsText={'BackgroundColor','yellow','FontSize',20};

E=epipolarBuildEFromRT(R1,T1,R2,T2);
[R12,T12]=computeRelativePoseFromRT(R1,T1,R2,T2);
fprintf('%s\n',titleText)
disp([E R12 T12])
draw3dcameraFromRT(R1,T1,optsDrawCamera{:})
hold on
draw3dcameraFromRT(R2,T2,optsDrawCamera{:})
plot3([T1(1) T2(1)],[T1(2) T2(2)],[T1(3) T2(3)])
text(T1(1),T1(2),T1(3),'1',optsText{:})
text(T2(1),T2(2),T2(3),'2',optsText{:})
hold off
axis square
axis equal
title(titleText)