clc
clear
close all
format short
addpath ./Functions

%Code that generate the Fig. 1 of 
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

K = @(s) (1-exp(-2*s))./(2*s); %kernel
phi = @(t) t.^(2.5).*exp(-t); %datum
d_phi = @(t) 2.5*t.^(2.5-1).*exp(-t) - t.^(2.5).*exp(-t); %derivative of the datum useful for the exact solution

Nt_up = 5; %minimum power of 2 for the time instants
Nt_down = 11; %maximum power of 2 for the time instants

alpha2 = 2; %grading of the mesh
err_2_Trap = zeros(1,Nt_up-Nt_down+1);
err_2_BDF2 = zeros(1,Nt_up-Nt_down+1);

for i = Nt_up:Nt_down

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha2; %time mesh

    g_app_Trap = backward_gcCQ_Trap(phi,K,t);
    g_app_BDF2 = backward_gcCQ_BDF2(phi,K,t);

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        k = 0 : floor(t(j+1)/2);
        g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
    end

    err_2_Trap(i-Nt_up+1) = norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf');
    err_2_BDF2(i-Nt_up+1) = norm(g_app_BDF2-g_ex,'inf')/norm(g_ex,'inf');
    
end

plot_error_gCQ({},err_2_Trap,2.^(-(Nt_up:Nt_down)),1);
plot_error_gCQ({},err_2_BDF2,2.^(-(Nt_up:Nt_down)),2);

alpha1 = 1; %grading of the mesh
err_1_Trap = zeros(1,Nt_up-Nt_down+1);
err_1_BDF2 = zeros(1,Nt_up-Nt_down+1);

for i = Nt_up:Nt_down

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha1;

    g_app_Trap = backward_gcCQ_Trap(phi,K,t);
    g_app_BDF2 = backward_gcCQ_BDF2(phi,K,t);

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        for k = 0 : floor(t(j+1)/2)
            g_ex(j) =  g_ex(j) + 2*d_phi(t(j+1)-2*k);
        end
    end

    err_1_Trap(i-Nt_up+1) = norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf');
    err_1_BDF2(i-Nt_up+1) = norm(g_app_BDF2-g_ex,'inf')/norm(g_ex,'inf');
end

plot_error_gCQ({'trap alpha=2','trap alpha=1'},err_1_Trap,2.^(-(Nt_up:Nt_down)),1);
plot_error_gCQ({'BDF2 alpha=2','BDF2 alpha=1'},err_1_BDF2,2.^(-(Nt_up:Nt_down)),2);