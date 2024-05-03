clc
clear
close all
format short
addpath ./Functions

%Code that generate the Fig. 5 (right)
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

p = 2.5;
K = @(s) (1-exp(-2*s))./(2*s); %kernel
phi = @(t) t.^p.*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^p.*exp(-t); %derivative of the datum useful for the exact solution

Nt_up = 5; %minimum power of 2 for the time instants
Nt_down = 7; %maximum power of 2 for the time instants

alpha2 = 2; %grading of the mesh

for i = Nt_up:Nt_down

    Nt = 2^i;
    t = ((0:Nt)/Nt).^alpha2; %time mesh

    %exact solution
    g_ex = zeros(Nt,1);
    for j = 1 : Nt
        k = 0 : floor(t(j+1)/2);
        g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
    end

    cont_j  = 0;
    set_Nq = 1:0.5:4;
    for j = set_Nq

        cont_j = cont_j + 1;
        Nq = floor(j*log(Nt)^2*Nt);

        g_app_Radau = backward_gcCQ_RK_with_Nq(phi,K,t,22,Nq);

        err_2_Radau(cont_j) = norm(g_app_Radau-g_ex,'inf')/norm(g_ex,'inf');

    end

    if i == 5
       
        semilogy(floor(set_Nq*Nt*log(Nt)^2),err_2_Radau,'or','LineWidth',1.7);
        hold on
        semilogy(floor(set_Nq*Nt*log(Nt)^2),err_2_Radau,'r','LineWidth',2);
        xline(floor(2*Nt*log(Nt)^2),'--','LineWidth',1.4)

    elseif i == 6

        semilogy(floor(set_Nq*Nt*log(Nt)^2),err_2_Radau,'*m','LineWidth',1.7);
        hold on
        semilogy(floor(set_Nq*Nt*log(Nt)^2),err_2_Radau,'m','LineWidth',2);
        xline(floor(2*Nt*log(Nt)^2),'--','LineWidth',1.4)

    elseif i == 7

        floor(set_Nq*Nt*log(Nt)^2)
        floor(2*Nt*log(Nt)^2)
        semilogy(floor(set_Nq*Nt*log(Nt)^2),err_2_Radau,'+b','LineWidth',1.7);
        hold on
        semilogy(floor(set_Nq*Nt*log(Nt)^2),err_2_Radau,'b','LineWidth',2);
        xline(floor(2*Nt*log(Nt)^2),'--','LineWidth',1.4)

    end

end