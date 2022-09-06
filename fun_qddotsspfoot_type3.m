function qddot = fun_qddotsspfoot_type3(x,u,dt)
%%%% SSP phase- 2   %%%%%%%%%%%%%%


global A B Nx Nu pert MI L m nx ny tx ty g  r lam  vars misc alp indic kc lamall xdata lamx lamy af acal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m4 = m(4);




MI4 = MI(4);


gamma61 = misc(5);
gamma62 = misc(6);
gamma71 = misc(7);
gamma72 = misc(8);
r4 = vars(1);
r7 = vars(2);
r4t = vars(3);
r4h = vars(4);



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


tht4 = x(1);
x1   = x(2);
y1   = x(3);
omg4 = x(4);
vhx = x(5);
vhy = x(6);
%{
alp4 = alp(4);
ax1 =  alp(8);
ay1 =  alp(9);
%}
%T1  = u(1);
T4  = u(1);

alp4 = alp(1);
ax1 =  alp(2);
ay1 =  alp(3);


x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding external forces
Ass1=[0,(-1),0,(-1);(-1),0,(-1),0;r4.*cos(tht4)+(-1).*r4t.*cos(gamma62+ ...
  tht4),(-1).*r4.*sin(tht4)+r4t.*sin(gamma62+tht4),r4.*cos(tht4)+( ...
  -1).*r4h.*cos(gamma61+tht4),(-1).*r4.*sin(tht4)+r4h.*sin(gamma61+ ...
  tht4);r4t.*cos(gamma62+tht4),0,r4h.*cos(gamma61+tht4),0]


bss1=[(-1).*m4.*(ax1+r4.*(omg4.^2.*cos(tht4)+alp4.*sin(tht4))),(-1).* ...
  m4.*(ay1+g+r4.*((-1).*alp4.*cos(tht4)+omg4.^2.*sin(tht4))),(-1).* ...
  alp4.*MI4+T4+m4.*r4.*((-1).*alp4.*r4+ay1.*cos(tht4)+g.*cos(tht4)+( ...
  -1).*ax1.*sin(tht4)),alp4.*MI4+m4.*((-1).*y1+r4.*sin(tht4)).*(ax1+ ...
  r4.*(omg4.^2.*cos(tht4)+alp4.*sin(tht4)))]



%bss1 = [bss11;bss12;bss13;bss14;bss15;bss16;bss17;bss18;bss19;bss20]

Assinv = inv(Ass1)
sol = Ass1\bss1'
%fvars = [sol ; det(Ass2)];



lam1 = sol(1);
lam2 = sol(2);
lam3 = sol(3);
lam4 = sol(4);

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


tht4 = x(1);
x1   = x(2);
y1   = x(3);
omg4 = x(4);
vhx  = x(5);
vhy  = x(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat =[m4.*r4.*sin(tht4),m4,0;(-1).*m4.*r4.*cos(tht4),0,m4;MI4+m4.* ...
  r4.^2,m4.*r4.*sin(tht4),(-1).*m4.*r4.*cos(tht4)]





%%%%%%%%%%%%%%%%%%%%%

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% phi is for finding qddot with ext lambda calculation 
phi=[lam2+lam4+(-1).*m4.*omg4.^2.*r4.*cos(tht4),lam1+lam3+(-1).*g.*m4+ ...
  (-1).*m4.*omg4.^2.*r4.*sin(tht4),T4+(-1).*(lam1+lam3).*r4.*cos( ...
  tht4)+g.*m4.*r4.*cos(tht4)+lam3.*r4h.*cos(gamma61+tht4)+lam1.* ...
  r4t.*cos(gamma62+tht4)+(lam2+lam4).*r4.*sin(tht4)+(-1).*lam4.* ...
  r4h.*sin(gamma61+tht4)+(-1).*lam2.*r4t.*sin(gamma62+tht4)]



%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NO contact

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For calculating lambda  
%{


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

%Gmat = Cn*invM*Cn';
%invG = inv(Gmat);

%nr  = Gd' - Cn*invM*rhs1';
%P1 = Gd' ;
%P2 = Cn*invM*rhs1';
%lam = Gmat\( Gd'- Cn*invM*rhs1')

lam = [lam1;lam2;lam3;lam4];
%(Cn'*lam + rhs1')
%qddot_invdy = invM*(Cn'*lam + rhs1');
%}

% qddot = qddot_invdy

phi;
qddot_invdy = invM*phi';
%alp = invM*phi'
%qddot_invdy(4) = 0
acal = qddot_invdy;
%af = qddot_invdy;
%qddot = [omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy;qddot_invdy]
 qddot = qddot_invdy
 %qddot = alp;
 lamx = lam2 + lam4;
 lamy = lam1 + lam3;
% lamx = lamt;
 %lamy = lamn;

end






 






