% date 5th Feb 2022
% Author quyinsong

% parameter
m11 = 25.8; m22 = 33.8; m33 = 2.76; m23 = 1.095; m32 = m23;
m0 = m22*m33-m23*m32;
syms u v r

c13 = -m22*v-m23*r; c31 = m22*v+m23*r;
c23 = m11*u; c32 = -m11*u;
d11 = -0.72253+1.32742*abs(u);
d22 = 0.88965+36.47287*abs(v)+0.805*abs(r);
d33 = 1.9-0.08*abs(v)+0.75*abs(r);
d23 = 7.25+0.845*abs(v)+3.45*abs(r);
d32 = -0.0313-3.95645*abs(v)-0.13*abs(r);
% cal fu
fu = (c13*r+d11*u)/m11;
% cal fv
fv =(-m23*c31*u-m23*c32*v+m33*c23*r+(m33*d22-m23*d32)*v+(m33*d23-m23*d33)*r)/m0;
% cal fr
fr = (m22*c31*u+m22*c32*v-m32*c23*r+(-m32*d22+m22*d32)*v+(-m32*d23+m22*d33)*r)/m0;
disp('fu')
vpa(fu)
disp('fv')
vpa(fv)
disp('fr')
vpa(fr)

