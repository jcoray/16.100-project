% MachSweepStudy.m
% Created by Team 2 - Jakob, Jose, and Austin on 11/29/18

L = 75000 * 9.81; % N
tc_p = 0.14;
S_ref = 127; % m^2
rho = 0.41351; % 10,000 m alt
b = 35.79; % m
e = 1;

% Part 1

M_inf = 0.78;
lambda = 25;

q_inf = 1/2 * rho * (300)^2 * M_inf^2;
cl_perp = L / (q_inf * S_ref * (cosd(lambda))^2);
Mdd = 0.86 - 0.75 * tc_p - 0.05 * cl_perp;

% Part 2

M_inf = 0.72541;

q_inf = 1/2 * rho * (300)^2 * M_inf^2
CL = L/(q_inf*S_ref);
p = [M_inf .75*tc_p-.86 0 .05*CL];
r = roots(p);
lamb = acosd(r)

D_i = (L/b)^2/(q_inf*pi*e)