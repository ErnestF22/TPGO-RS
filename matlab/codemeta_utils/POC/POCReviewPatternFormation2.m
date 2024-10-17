function POCReviewPatternFormation2
%
%Check of CDC2018 submission, eq. (17).
syms kv gp kp ci di mu
h=@(ci,di) mu^2-kv*mu*(ci+i*di)-(-gp+kp*(ci+i*di))
collect(expand(h(ci,di)*h(ci,-di)),mu)
simplify(ci^2*kp^2 - 2*ci*gp*kp + di^2*kp^2 + gp^2 - ((gp-kp*ci)^2+kp^2*di^2))
%might need debugging