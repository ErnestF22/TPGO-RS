function dx=bearingControlField(y,yg,funsTheta)
[v,theta]=cnormalize(sphere_log(y,yg));
dx=-funsTheta.f(theta)*v;
