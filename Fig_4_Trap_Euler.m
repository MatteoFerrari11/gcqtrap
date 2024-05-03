clc
clear
close all
format short
addpath ./Functions

%Code that generate the Fig. 4 of 
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

%%% p = 1.5 %%%%

p = 1.5;
K = @(z) (1-exp(-2*z))./(2*z);
phi = @(t) t.^(p).*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^(p).*exp(-t); %derivative of the datum useful for the exact solution

Nt_up = 5; %minimum power of 2 for the time instants
Nt_down = 11; %maximum power of 2 for the time instants

alpha2 = 2; %grading of the mesh
err_2_Euler = zeros(1,Nt_up-Nt_down+1);
err_2_Trap = zeros(1,Nt_up-Nt_down+1);

for i = Nt_up:Nt_down

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha2; %time mesh

    g_app_Trap = backward_gcCQ_Trap(phi,K,t);
    g_app_Euler = backward_gcCQ_Euler(phi,K,t);

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        k = 0 : floor(t(j+1)/2);
        g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
    end

    err_2_Trap(i-Nt_up+1) = norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf');
    err_2_Euler(i-Nt_up+1) = norm(g_app_Euler-g_ex,'inf')/norm(g_ex,'inf');
end

plot_error_gCQ({},err_2_Trap,2.^(-(Nt_up:Nt_down)),1);
hold on
plot_error_gCQ({'trap alpha=2','Euler alpha=2'},err_2_Euler,2.^(-(Nt_up:Nt_down)),1);
title('p=1.5')

%%% p = 2.5 %%%%

p = 2.5;
K = @(z) (1-exp(-2*z))./(2*z);
phi = @(t) t.^(p).*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^(p).*exp(-t); %derivative of the datum useful for the exact solution

Nt_up = 5; %minimum power of 2 for the time instants
Nt_down = 11; %maximum power of 2 for the time instants

alpha2 = 2; %grading of the mesh
err_2_Euler = zeros(1,Nt_up-Nt_down+1);
err_2_Trap = zeros(1,Nt_up-Nt_down+1);

for i = Nt_up:Nt_down

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha2; %time mesh

    g_app_Trap = backward_gcCQ_Trap(phi,K,t);
    g_app_Euler = backward_gcCQ_Euler(phi,K,t);

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        k = 0 : floor(t(j+1)/2);
        g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
    end

    err_2_Trap(i-Nt_up+1) = norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf');
    err_2_Euler(i-Nt_up+1) = norm(g_app_Euler-g_ex,'inf')/norm(g_ex,'inf');
end

plot_error_gCQ({},err_2_Trap,2.^(-(Nt_up:Nt_down)),2);
hold on
plot_error_gCQ({'trap alpha=2','Euler alpha=2'},err_2_Euler,2.^(-(Nt_up:Nt_down)),2);
title('p=2.5')

