% Script to load the results into Octave
clear all
RHO_all = load("../RUN/rho.out");
U_all = load("../RUN/u.out");
P_all = load("../RUN/P.out");
gam_all = load("../RUN/gamma.out");
Pc_all = load("../RUN/Pc.out");
E_all = load("../RUN/E.out");
A_all = load("../RUN/A.out");
%

% #######################################################
% Read some of the parameters from init.params
nr = 200
outf = 200
% #######################################################
total = size(RHO_all)
nt = total(1)/nr
%
%
for j=1:nt
	for i=1:nr
		r(i,j) = RHO_all(nr*(j-1) + i,1);
		rho(i,j) = RHO_all(nr*(j-1) + i,2);
		u(i,j) = U_all(nr*(j-1) + i,2);
		P(i,j) = P_all(nr*(j-1) + i,2);
		gam(i,j) = gam_all(nr*(j-1) + i,2);
		Pc(i,j) = Pc_all(nr*(j-1) + i,2);
		E(i,j) = E_all(nr*(j-1) + i,2);
		alpha(i,j) = A_all(nr*(j-1) + i,2);
%
	endfor
endfor
little_e = E./rho .- 0.5.*u.**2;

Pinf = P(nr,1);
PovPinf = (1/Pinf).*P;
rhonorm = rho./max(rho(1,1),rho(nr,1));
unorm = u./max(max(u(nr,:)),1);
%PovPinfi = (1/Pinf).*Pi;
%
T_all = load("../RUN/interface.out");
Pb_all = load("../RUN/pint.out");
T = T_all(:,1);
Rb = T_all(:,2);
Pbubble = Pb_all(:,2);
for i=1:nt
	t(i) = T_all(outf*(i),1);
	R(i) = T_all(outf*(i),2);
endfor
