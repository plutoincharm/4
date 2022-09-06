function qddot = fun_qddotssp_leganksp(x,u,dt)
%%%% SSP phase-1 of stance leg , leg and ankle   %%%%%%%%%%%%%%


global A B Nx Nu pert MI L m nx ny tx ty g  r lam  vars misc alp indic kc lamall xdata lamx lamy af acal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1 = m(1);
m2 = m(2);
m3 = m(3);
m4 = m(4);


L2 = L(1);
L3 = L(2);



r2 = r(1);
r3 = r(2);


MI1 = MI(1);
MI2 = MI(2);
MI3 = MI(3);
MI4 = MI(4);

R2 = misc(1);
R3 = misc(2);
%R2 = 0.4028;
%R3 = 0.4061;
del2 = misc(3);
del3 = misc(4);
%del2 = deg2rad(1.142);
%del3 = deg2rad(1.20);
gamma61 = misc(5);
gamma62 = misc(6);
gamma71 = misc(7);
gamma72 = misc(8);
r4 = vars(1);
r7 = vars(2);
r4t = vars(3);
r4h = vars(4);

%M= 0;
d1 = 0.1081;
d3 = 0.0682;
h4 = 0.0437;
h2 = 0.0722;
r4c =  vars(5)   %0.1638;
% 0.0872 new value considering lycop = 0
gamma43 =  vars(6) % 2.3752 ;  
% 2.6481 new value considering lycop = 0

%{
tht4 = xdata(1);
x1   = xdata(2);
y1   = xdata(3);
omg4 = xdata(4);
vhx = xdata(5);
vhy = xdata(5);
alp4 = alp(4);
ax1 =  alp(8);
ay1 =  alp(9);
%T1  = u(1);
T4  = u(1);
%}


tht3 = x(1);
tht4 = x(2);
x1   = x(3);
y1   = x(4);
omg3 = x(5);
omg4 = x(6);
vhx = x(7);
vhy = x(8);

%{
alp4 = alp(4);
ax1 =  alp(8);
ay1 =  alp(9);
%}
%T2  = u(1);
T3  = u(1);
T4  = u(2);




alp3 = alp(1);
alp4 = alp(2);
ax1 =  alp(3);
ay1 =  alp(4);
lam1 = 44;
lam2 = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding external forces
%{
Ass1=[0,(-1),0,(-1);(-1),0,(-1),0;r4.*cos(tht4)+(-1).*r4t.*cos(gamma62+ ...
  tht4),(-1).*r4.*sin(tht4)+r4t.*sin(gamma62+tht4),r4.*cos(tht4)+( ...
  -1).*r4h.*cos(gamma61+tht4),(-1).*r4.*sin(tht4)+r4h.*sin(gamma61+ ...
  tht4);(m3+m4).^(-1).*(m4.*r4t.*cos(gamma62+tht4)+m3.*((-1).*L3.* ...
  cos(tht3)+r3.*cos(tht3)+(-1).*r4.*cos(tht4)+r4t.*cos(gamma62+tht4) ...
  )),0,(m3+m4).^(-1).*(m4.*r4h.*cos(gamma61+tht4)+m3.*((-1).*L3.* ...
  cos(tht3)+r3.*cos(tht3)+(-1).*r4.*cos(tht4)+r4h.*cos(gamma61+tht4) ...
  )),0]

bss1=[(-1).*m3.*(ax1+r3.*(omg3.^2.*cos(tht3)+alp3.*sin(tht3)))+(-1).* ...
  m4.*(ax1+L3.*(omg3.^2.*cos(tht3)+alp3.*sin(tht3))+r4.*(omg4.^2.* ...
  cos(tht4)+alp4.*sin(tht4))),(-1).*m3.*(ay1+g+r3.*((-1).*alp3.*cos( ...
  tht3)+omg3.^2.*sin(tht3)))+(-1).*m4.*(ay1+g+(-1).*alp4.*r4.*cos( ...
  tht4)+L3.*((-1).*alp3.*cos(tht3)+omg3.^2.*sin(tht3))+omg4.^2.*r4.* ...
  sin(tht4)),(-1).*alp4.*MI4+T4+m4.*r4.*((-1).*alp4.*r4+(-1).*alp3.* ...
  L3.*cos(tht3+(-1).*tht4)+ay1.*cos(tht4)+g.*cos(tht4)+L3.*omg3.^2.* ...
  sin(tht3+(-1).*tht4)+(-1).*ax1.*sin(tht4)),alp3.*MI3+alp4.*MI4+ ...
  m3.*((-1).*y1+r3.*sin(tht3)).*(ax1+r3.*(omg3.^2.*cos(tht3)+alp3.* ...
  sin(tht3)))+m3.*m4.*(m3+m4).^(-1).*(L3.*cos(tht3)+(-1).*r3.*cos( ...
  tht3)+r4.*cos(tht4)).*(ay1+g+r3.*((-1).*alp3.*cos(tht3)+omg3.^2.* ...
  sin(tht3)))+(-1).*m3.*m4.*(m3+m4).^(-1).*(L3.*cos(tht3)+(-1).*r3.* ...
  cos(tht3)+r4.*cos(tht4)).*(ay1+g+(-1).*alp4.*r4.*cos(tht4)+L3.*(( ...
  -1).*alp3.*cos(tht3)+omg3.^2.*sin(tht3))+omg4.^2.*r4.*sin(tht4))+ ...
  m4.*((-1).*y1+L3.*sin(tht3)+r4.*sin(tht4)).*(ax1+L3.*(omg3.^2.* ...
  cos(tht3)+alp3.*sin(tht3))+r4.*(omg4.^2.*cos(tht4)+alp4.*sin(tht4) ...
  ))]



%bss1 = [bss11;bss12;bss13;bss14;bss15;bss16;bss17;bss18;bss19;bss20]

Assinv = inv(Ass1)
sol = Ass1\bss1'
%fvars = [sol ; det(Ass2)];



lam1 = sol(1);
lam2 = sol(2);
lam3 = sol(3);
lam4 = sol(4);
%}
%{
lamn = sol(1) + sol(3);
lamt = sol(2) + sol(4);
lam1 = lamn/2;
lam3 = lamn/2;
lam2 = lamt/2;
lam4 = lamt/2;
%}



%{
if (indic ==15)
    lam1 = lamall(1,kc);
    lam2 = lamall(2,kc);
    lam3 = lamall(3,kc);  
    lam4 = lamall(4,kc);
    
end
%}
%taking torque values from calculation
%{
T2  = sol(5);
T3  = sol(6);
T4  = sol(7);
T5  = sol(8);
T6  = sol(9);
T7  = sol(10);
%}  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat=[L3.*m4.*sin(tht3)+m3.*r3.*sin(tht3),m4.*r4.*sin(tht4),m3+m4,0;( ...
  -1).*L3.*m4.*cos(tht3)+(-1).*m3.*r3.*cos(tht3),(-1).*m4.*r4.*cos( ...
  tht4),0,m3+m4;L3.^2.*m4+MI3+m3.*r3.^2,L3.*m4.*r4.*cos(tht3+(-1).* ...
  tht4),L3.*m4.*sin(tht3)+m3.*r3.*sin(tht3),(-1).*L3.*m4.*cos(tht3)+ ...
  (-1).*m3.*r3.*cos(tht3);L3.*m4.*r4.*cos(tht3+(-1).*tht4),MI4+m4.* ...
  r4.^2,m4.*r4.*sin(tht4),(-1).*m4.*r4.*cos(tht4)]





%%%%%%%%%%%%%%%%%%%%%

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% phi is for finding qddot with ext lambda calculation 

phi=[lam2+(-1).*L3.*m4.*omg3.^2.*cos(tht3)+(-1).*m3.*omg3.^2.*r3.*cos( ...
  tht3)+(-1).*m4.*omg4.^2.*r4.*cos(tht4),lam1+(-1).*g.*m3+(-1).*g.* ...
  m4+(-1).*L3.*m4.*omg3.^2.*sin(tht3)+(-1).*m3.*omg3.^2.*r3.*sin( ...
  tht3)+(-1).*m4.*omg4.^2.*r4.*sin(tht4),T3+(-1).*T4+(-1).*L3.* ...
  lam1.*cos(tht3)+g.*L3.*m4.*cos(tht3)+g.*m3.*r3.*cos(tht3)+L3.* ...
  lam2.*sin(tht3)+(-1).*L3.*m4.*omg4.^2.*r4.*sin(tht3+(-1).*tht4), ...
  T4+(-1).*lam1.*r4.*cos(tht4)+g.*m4.*r4.*cos(tht4)+lam1.*r4c.*cos( ...
  gamma43+tht4)+L3.*m4.*omg3.^2.*r4.*sin(tht3+(-1).*tht4)+lam2.*r4.* ...
  sin(tht4)+(-1).*lam2.*r4c.*sin(gamma43+tht4)]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NO contact


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For calculating lambda  
%{
Cn=[0,(-1),0,(-1);(-1),0,(-1),0;d1+L3.*cos(tht3),h2+(-1).*L3.*sin( ...
  tht3),(-1).*d3+L3.*cos(tht3),h4+(-1).*L3.*sin(tht3);(-1).*d1+r4.* ...
  cos(tht4)+(-1).*r4t.*cos(gamma62+tht4),(-1).*h2+(-1).*r4.*sin( ...
  tht4)+r4t.*sin(gamma62+tht4),d3+r4.*cos(tht4)+(-1).*r4h.*cos( ...
  gamma61+tht4),(-1).*h4+(-1).*r4.*sin(tht4)+r4h.*sin(gamma61+tht4) ...
  ]

Gd=[(-1).*L3.*omg3.^2.*sin(tht3)+(-1).*omg4.^2.*r4.*sin(tht4)+ ...
  omg4.^2.*r4t.*sin(gamma62+tht4),(-1).*L3.*omg3.^2.*cos(tht3)+(-1) ...
  .*omg4.^2.*r4.*cos(tht4)+omg4.^2.*r4t.*cos(gamma62+tht4),(-1).* ...
  L3.*omg3.^2.*sin(tht3)+(-1).*omg4.^2.*r4.*sin(tht4)+omg4.^2.*r4h.* ...
  sin(gamma61+tht4),(-1).*L3.*omg3.^2.*cos(tht3)+(-1).*omg4.^2.*r4.* ...
  cos(tht4)+omg4.^2.*r4h.*cos(gamma61+tht4)]

rhs1=[(-1).*L3.*m4.*omg3.^2.*cos(tht3)+(-1).*m3.*omg3.^2.*r3.*cos(tht3)+...
  (-1).*m4.*omg4.^2.*r4.*cos(tht4),(-1).*g.*m3+(-1).*g.*m4+(-1).* ...
  L3.*m4.*omg3.^2.*sin(tht3)+(-1).*m3.*omg3.^2.*r3.*sin(tht3)+(-1).* ...
  m4.*omg4.^2.*r4.*sin(tht4),T3+(-1).*T4+g.*L3.*m4.*cos(tht3)+g.* ...
  m3.*r3.*cos(tht3)+(-1).*L3.*m4.*omg4.^2.*r4.*sin(tht3+(-1).*tht4), ...
  T4+g.*m4.*r4.*cos(tht4)+L3.*m4.*omg3.^2.*r4.*sin(tht3+(-1).*tht4) ...
  ]




%}

%{
%%%%%%%%%%%%%%%%%%%%%%%
inv(Mmat);
qddot = Mmat\phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhs1 = rhs;
%Cn  = Cmat;
%rhs2 = Gd';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%rhs2 -gn
%%rhs1 - phi
%Mmat
invM = inv(Mmat)

%Gmat = Cn*invM*Cn'
%invG = inv(Gmat)

%nr  = Gd' - Cn*invM*rhs1';
%P1 = Gd' ;
%P2 = Cn*invM*rhs1';
%lam = Gmat\( Gd'- Cn*invM*rhs1')

%lam = [lam1;lam2;lam3;lam4];
%(Cn'*lam + rhs1')
%qddot_invdy = invM*(Cn'*lam + rhs1');


% qddot = qddot_invdy

%phi;
qddot_invdy = invM*phi';
%acal = qddot_invdy
%alp = invM*phi'
%qddot_invdy(4) = 0
af = qddot_invdy;
%qddot = [omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy;qddot_invdy]
 qddot = qddot_invdy
 %qddot = alp;
% lamx = lam(2) + lam(4);
% lamy = lam(1) + lam(3);
% lamx = lam2 + lam4;
%lamy = lam1 + lam3;
% lamx = lamt;
 %lamy = lamn;

end






 






