clc
clear
close all
format short
addpath ./Functions

%Code that generates Table 1 at pag. 23 of
% [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule]

p = 2.5;
K = @(z) (1-exp(-2*z))./(2*z);
phi = @(t) t.^(p).*exp(-t); %datum
d_phi = @(t) p*t.^(p-1).*exp(-t) - t.^(p).*exp(-t); %derivative of the datum useful for the exact solution

alpha2 = 2; %grading of the mesh
time_Euler = zeros(1,4);
time_Trap = zeros(1,4);

cont = 1;
for err = [10^(-1) 10^(-2) 10^(-3) 10^(-4) 10^(-5)]  

    for i = 1:20
        
        tic
        Nt = 2^i;
        t = ((0:Nt)/Nt).^alpha2; %time mesh

        g_app_Trap = backward_gcCQ_Trap(phi,K,t);

        %exact solution
        g_ex = zeros(Nt,1);
        for j = 1 : Nt
            k = 0 : floor(t(j+1)/2);
            g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
        end
        time_Trap_cor = toc;
        if norm(g_app_Trap-g_ex,'inf')/norm(g_ex,'inf') < err
            time_Trap(cont) = time_Trap_cor;
            break
        end
    end

    for i = 1:20

        tic
        Nt = 2^i;
        t = ((0:Nt)/Nt).^alpha2; %time mesh

        g_app_Euler = backward_gcCQ_Euler(phi,K,t);

        %exact solution
        g_ex = zeros(Nt,1);
        for j = 1 : Nt
            k = 0 : floor(t(j+1)/2);
            g_ex(j) =  g_ex(j) + sum(2*d_phi(t(j+1)-2*k));
        end
        time_Euler_cor = toc; 

        if norm(g_app_Euler-g_ex,'inf')/norm(g_ex,'inf') < err
            time_Euler(cont) = time_Euler_cor;
            break
        end
    end
     cont = cont+1;
 end

time_Trap
time_Euler
