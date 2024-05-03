clc
clear
close all
format short
addpath ./Functions

%Code that generate the Fig. 5 (left)
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

p = 2.5;
K = @(s) (1-exp(-2*s))./(2*s); %kernel
phi = @(t) t.^p.*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^p.*exp(-t); %derivative of the datum useful for the exact solution

Nt_up = 3; %minimum power of 2 for the time instants
Nt_down = 8; %maximum power of 2 for the time instants

alpha2 = 2; %grading of the mesh
err_2_Radau = zeros(1,Nt_up-Nt_down+1);
err_2_Trap = zeros(1,Nt_up-Nt_down+1);

time_Trap = 0;
time_Radau = 0;

for i = Nt_up:Nt_down

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha2; %time mesh

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        k = 0 : floor(t(j+1)/2);
        g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
    end

    g_app_Radau = backward_gcCQ_RK(phi,K,t,22);
    
    err_2_Radau(i-Nt_up+1) = norm(g_app_Radau-g_ex,'inf')/norm(g_ex,'inf');

    
    g_app_Trap = backward_gcCQ_Trap(phi,K,t);
    
    err_2_Trap(i-Nt_up+1) = norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf');
end

plot_error_gCQ({},err_2_Radau,2.^(-(Nt_up:Nt_down)),1);
hold on
plot_error_gCQ({'Radau IIA alpha=2','Trap alpha=2'},err_2_Trap,2.^(-(Nt_up:Nt_down)),1);