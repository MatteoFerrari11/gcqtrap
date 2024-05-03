clc
clear
close all
format long
addpath ./Functions

%Code that generates Table 2 at pag. 29 of
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

p = 2.5;
K = @(z) (1-exp(-2*z))./(2*z);
phi = @(t) t.^(p).*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^(p).*exp(-t); %derivative of the datum useful for the exact solution

alpha = 2; %grading of the mesh
err_Trap = zeros(1,2);
err_RK = zeros(1,2);

time_Trap = zeros(1,2);
time_RK = zeros(1,2);

cont = 1;

for i = [7 9]

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha; %time mesh

    tic
    g_app_Trap = backward_gcCQ_Trap(phi,K,t);
    time_Trap(cont) = toc;

    Nq = 2*floor(Nt*log(Nt)^2);

    tic
    g_app_RK = backward_gcCQ_RK(phi,K,t,22);
    time_RK(cont) = toc;

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        k = 0 : floor(t(j+1)/2);
        g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
    end

    err_Trap(cont) = norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf');
    err_RK(cont) = norm(g_app_RK-g_ex,'inf')/norm(g_ex,'inf');
    cont = cont+1;
end

err_Trap
err_RK

time_Trap
time_RK

