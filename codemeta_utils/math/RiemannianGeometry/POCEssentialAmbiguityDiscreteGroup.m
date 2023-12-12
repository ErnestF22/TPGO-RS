function POCEssentialAmbiguityDiscreteGroup
%syms theta
theta=randn;
Rzt=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
Hzt={cat(3,Rzt,Rzt)};

I=eye(3);
Rxpi=diag([1;-1;-1]);
Rypi=diag([-1;1;-1]);
Rzpi=diag([-1;-1;1]);
Hz={cat(3,eye(3),eye(3)), cat(3,eye(3),Rzpi)};
Hx={cat(3,eye(3),eye(3)), cat(3,Rxpi,Rxpi)};
Hzx=cartesianGroupProduct(Hz,Hx);
Hzx2=cartesianGroupProduct(Hzx,Hzx);
HztHzx=cartesianGroupProduct(Hzt,Hzx);
HzxHzt=cartesianGroupProduct(Hzx,Hzt);


disp('Hzx')
groupDisplay(Hzx)
disp('cellUnique(Hzx2)')
groupDisplay(cellUnique(Hzx2))
disp('[I Rxpi I Rxpi; I Rxpi Rzpi Rypi]')
disp([I Rxpi I Rxpi; I Rxpi Rzpi Rypi])
disp('HztHzx')
groupDisplay(HztHzx)
disp('HzxHzt')
groupDisplay(HzxHzt)

disp('Cayley Table')
cayleyTable(Hzx)

function H3=cartesianGroupProduct(H1,H2)
NH1=length(H1);
NH2=length(H2);
H3=cell(1,NH1*NH2);

cnt=0;
for iH1=1:NH1
    for iH2=1:NH2
        cnt=cnt+1;
        H3{cnt}=multRotArray(H1{iH1},H2{iH2});
    end
end

function H=cellUnique(H)
NH=length(H);
flagRemove=false(1,NH);
for iH=1:NH
    for jH=1:iH-1
        if all(H{iH}(:)==H{jH}(:))
            flagRemove(iH)=true;
            break
        end
    end
end
H=H(~flagRemove);

function groupDisplay(H)
c1=cellfun(@(x) x(:,:,1),H,'UniformOutput',false);
c2=cellfun(@(x) x(:,:,2),H,'UniformOutput',false);
disp([c1{:};c2{:}])

function cayleyTable(H)
NH=length(H);
for iH=1:NH
    for jH=1:NH
        Hprod=multRotArray(H{iH},H{jH});
        idxHprod=find(cellfun(@(x) all(x(:)==Hprod(:)), H),1,'first');
        fprintf( '%d ',idxHprod)
    end
    fprintf('\n')
end

