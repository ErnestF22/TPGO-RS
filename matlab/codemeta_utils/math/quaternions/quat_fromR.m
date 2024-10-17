function q=quat_fromR(R)
[r,theta]=logrot(R);
q=[cos(theta/2); r.*([1;1;1]*sin(theta/2))];
