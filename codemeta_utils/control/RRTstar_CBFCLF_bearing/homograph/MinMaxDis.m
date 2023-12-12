function [sMin,sMax,zMax,zMin] = MinMaxDis(y,L,flag_test)
sMin = inf*ones(size(L,2),1);
sMax = zeros(size(L,2),1);
zMax = sMax;
zMin = sMin;
for i=1:size(L,2)
    for j=1:size(y,2)
        t = 1/norm(L(:,i)-y(:,j));
        if t<=sMin(i)
            if ~flag_test
                sMin(i)=t;
                zMax(i) = 1/t;
            else
                sMin(i)= 1;
            end
        end
        if t>sMax(i)
            if ~flag_test
                sMax(i)=t;
                zMin(i) = 1/t;
            else
                sMax(i)= 1.35;
            end
        end
    end
end

end