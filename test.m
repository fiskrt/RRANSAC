

%syms f(h)
%assume(h >= 0.001)
eps = 0.9;
eps_i = 0.7;
del_i = 0.1;

%eqn = eps*(eps_i/del_i)^h + (1-eps)*((1-eps_i)/(1-del_i))^h == 1;
% solve(eqn,h)
% 
syms h
 eq = eps*(del_i/eps_i).^h + (1-eps)*((1-del_i)/(1-eps_i)).^h;
% 
% 
sol = vpasolve(eq==1, h, [0.001, 1500])

fun = @(h) eps*(del_i/eps_i).^h + (1-eps)*((1-del_i)/(1-eps_i)).^h;
h = linspace(-2, 1.1);

figure
plot(h, fun(h))
