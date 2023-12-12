function dx=cartModel(x,u)
dx=[cos(x(3)) 0; sin(x(3)) 0; 0 1]*u;
