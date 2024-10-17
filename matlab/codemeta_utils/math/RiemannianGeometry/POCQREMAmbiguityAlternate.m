function POCQREMAmbiguityAlternate

disp('Constraints for zero coordinates')
Pz=diag([1 1 0]);
Si=sym('S',[3 3]);
disp(Si*Pz==Pz*Si)

disp('Constraints on angles')
e3hat=[0 -1 0; 1 0 0; 0 0 0];
Rxpi=diag([1;-1;-1]);
Rypi=diag([-1;1;-1]);
Rzpi=diag([-1;-1;1]);
syms theta theta1 theta2
Rz=@(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

h=getSeqCombSet('+-', 3);
for iseq=1:8
    s=h();
    fprintf('Case "%s"\n\n',s);
    if s(1)=='+'
        sign=+1;
    else
        sign=-1;
    end
    switch s
        case '+++'
            S1=@(theta1) Rz(theta1);
            S2=@(theta2) Rz(theta2);
            S1Str='S1Solve=Rz(theta);';
            S2Str='S2Solve=Rz(theta);';
        case '+-+'
            S1=@(theta1) Rypi*Rz(theta1);
            S2=@(theta2) Rz(theta2);
            S1Str='S1Solve=[];';
            S2Str='S2Solve=[];';
        case '++-'
            S1=@(theta1) Rz(theta1);
            S2=@(theta2) Rypi*Rz(theta2);
            S1Str='S1Solve=[];';
            S2Str='S2Solve=[];';
        case '+--'
            S1=@(theta1) Rypi*Rz(theta1);
            S2=@(theta2) Rypi*Rz(theta2);
            S1Str='S1Solve=Rypi*Rz(theta);';
            S2Str='S2Solve=Rypi*Rzpi*Rz(theta);';
        case '--+'
            S1=@(theta1) Rypi*Rz(theta1);
            S2=@(theta2) Rz(theta2);
            S1Str='S1Solve=[];';
            S2Str='S2Solve=[];';
        case '-+-'
            S1=@(theta1) Rz(theta1);
            S2=@(theta2) Rypi*Rz(theta2);
            S1Str='S1Solve=[];';
            S2Str='S2Solve=[];';
        case '-++'
            S1=@(theta1) Rz(theta1);
            S2=@(theta2) Rz(theta2);
            S1Str='S1Solve=Rz(theta);';
            S2Str='S2Solve=Rzpi*Rz(theta);';
        case '---'
            S1=@(theta1) Rypi*Rz(theta1);
            S2=@(theta2) Rypi*Rz(theta2);
            S1Str='S1Solve=Rypi*Rz(theta);';
            S2Str='S2Solve=Rypi*Rz(theta);';
    end
    
    eval(S1Str);
    eval(S2Str);
    
    fprintf('\tS1(theta)\n\n')
    S=S1(theta);
    disp(S)
    if (s(2)=='+' && S(3,3)~=1) || (s(2)=='-' && S(3,3)~=-1)
        error('Wrong s_i')
    end
    fprintf('\tCheck det(S1(theta))=1\n\n')
    d=simplify(det(S1(theta)));
    disp(d)
    if d~=1
        error('Wrong determinant')
    end
    
    fprintf('\tS2(theta)\n\n')
    S=S2(theta);
    disp(S)
    if (s(3)=='+' && S(3,3)~=1) || (s(3)=='-' && S(3,3)~=-1)
        error('Wrong s_i')
    end
    fprintf('\tCheck det(S1(theta))=1\n\n')
    d=simplify(det(S2(theta)));
    disp(d)
    if d~=1
        error('Wrong determinant')
    end

    c1=S1(theta1)*e3hat;
    c2=sign*e3hat*S2(theta2);
    fprintf('\tConstraints\n\n')
    disp([c1([1 2 4 5]).'==c2([1 2 4 5]).'])

    if ~isempty(S1Solve)
        fprintf('\tSolutions\n%s\n%s\n',S1Str,S2Str)
        fprintf('\tCheck S1Solve*e3hat-e3hat*S2Solve=0\n\n')
        err=S1Solve*e3hat-sign*e3hat*S2Solve;
        errEq=err==0;
        disp(err)
        if ~all(errEq(:))
            error('Constraints not satisfied')
        end
    else
        fprintf('\tCheck: case not possible\n\n')
        %pause
    end
end        
