function R = test_cumulative
figure(3)
A1 = [1 0;0 1];
b1 = -2*[1;1];

A2 = [4 1;3 1];
b2 = -5*[1;1];


A3 = [1 1];
b3 = 1;

x = -10:0.05:10;

nbColors=32;
colors=parula(nbColors);
gs1=getSeqCombSet([0 1],2);
cnt=0;
for is1=1:4
    z1=gs1()';
    gs2=getSeqCombSet([0 1],2);
    for is2=1:4
        z2=gs2()';
        gs3=getSeqCombSet([0 1],1);
        for is3=1:2
            z3=gs3()';
            cnt = cnt+1;
             
%             [A,b]= cumulativeCoeff(A1,b1,A2,b2,A3,b3,z1,z2,z3);
%             R(cnt).A = A;
%             R(cnt).b = b;
%             R(cnt).z1 = z1;
%             R(cnt).z2 = z2;
%             R(cnt).z3 = z3;
          
            Y = funImage(x,x,@(x) network(A1,b1,A2,b2,A3,b3,z1,z2,z3,x),'method','surf','methodOpts',{'FaceColor',colors(mod(cnt-1,nbColors)+1,:)});
            if ~all(all(isnan(Y)))
                legendInfo{cnt} = ['Z = ' num2str([z1',z2',z3])];
            else
                legendInfo{cnt} = [''];
            end
            hold on
            
        end
    end
end
legend(legendInfo)
% hold on
% plot3(2,7,1,'c*')
% hold on
% plot3(7,2,1,'m*')
% hold on

end

