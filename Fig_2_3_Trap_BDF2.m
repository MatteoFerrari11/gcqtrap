clc
clear
close all
format short
addpath ./Functions

%Code that generate the Fig. 2 and Fig. 3 of
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

p = 5/2;
K = @(s) (1-exp(-2*s))./(2*s); %kernel
phi = @(t) t.^p.*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^p.*exp(-t); %derivative of the datum useful for the exact solution

Nt = 64;

%%%%%%%% alpha = 1 %%%%%%%%
alpha1 = 1; %grading of the mesh

t_1 = ((0:Nt)/Nt).^alpha1;

g_app_Trap_1 = backward_gcCQ_Trap(phi,K,t_1);
g_app_BDF2_1 = backward_gcCQ_BDF2(phi,K,t_1);

%exact solution
g_ex_1 = zeros(Nt,1);
for j = 1 : Nt
    for k = 0 : floor(t_1(j+1)/2)
        g_ex_1(j) =  g_ex_1(j) + 2*d_phi(t_1(j+1)-2*k);
    end
end

figure(1)
plot(t_1,abs([0 ; g_ex_1]-[0 ; g_app_Trap_1]),'-o','linewidth',1.5,'Markersize',6)
xlabel('t','Fontsize',18)
ylabel('Absolute Error','Fontsize',18)
title('Trapezoidal rule alpha=1')
xtick = [0, 0.2, 0.4, 0.6, 0.8, 1];
xticks(xtick)
ax = gca;
ax.FontSize = 14;
ytick = [0, 0.04, 0.08, 0.12, 0.16, 0.20]*10^(-2);
yticks(ytick)

figure(2)
plot(t_1,abs([0 ; g_ex_1]-[0 ; g_app_BDF2_1]),'-o','linewidth',1.5,'Markersize',6)
xlabel('t','Fontsize',18)
ylabel('Absolute Error','Fontsize',18)
title('BDF2 alpha=1')
xtick = [0, 0.2, 0.4, 0.6, 0.8, 1];
xticks(xtick)
ax = gca;
ax.FontSize = 14;
ytick = [0, 0.1, 0.2, 0.3, 0.4, 0.5]*10^(-2);
yticks(ytick)

%%%%%%%% alpha = 2 %%%%%%%%
alpha2 = 2; %grading of the mesh

t_2 = ((0:Nt)/Nt).^alpha2; %time mesh

g_app_Trap_2 = backward_gcCQ_Trap(phi,K,t_2);
g_app_BDF2_2 = backward_gcCQ_BDF2(phi,K,t_2);

%exact solution
g_ex_2 = zeros(Nt,1);

for j = 1 : Nt
    k = 0 : floor(t_2(j+1)/2);
    g_ex_2(j) =  g_ex_2(j) + sum(2*d_phi(t_2(j+1)-2*k));
end

figure(3)
plot(t_2,abs([0 ; g_ex_2]-[0 ; g_app_Trap_2]),'-o','linewidth',1.5,'Markersize',6)
xlabel('t','Fontsize',18)
ylabel('Absolute Error','Fontsize',18)
title('Trapezoidal rule alpha=2')
xtick = [0, 0.2, 0.4, 0.6, 0.8, 1];
xticks(xtick)
ax = gca;
ax.FontSize = 14;
ytick = [0, 0.4, 0.8, 1.2, 1.6, 2]*10^(-4);
yticks(ytick)

figure(4)
plot(t_2,abs([0 ; g_ex_2]-[0 ; g_app_BDF2_2]),'-o','linewidth',1.5,'Markersize',6)
xlabel('t','Fontsize',18)
ylabel('Absolute Error','Fontsize',18)
title('BDF2 alpha=2')
xtick = [0, 0.2, 0.4, 0.6, 0.8, 1];
xticks(xtick)
ax = gca;
ax.FontSize = 14;
ytick = [0, 0.4, 0.8, 1.2, 1.6, 2]*10^(-3);
yticks(ytick)

