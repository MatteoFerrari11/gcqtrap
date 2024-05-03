function g_n = backward_gcCQ_Euler(phi,K_op,t)

% function that compute the backward generalised convolution quadrature based on the backward Euler method 
% (for details see [Lopez Fernandez, Sauter, Generalized convolution quadrature with variable time stepping.
% Part II: Algorithm and numerical results]) to solve convolution equation k * g = phi

%input:
%      - phi : given function datum 
%      - K_op : the laplace transform of the kernel k
%      - t : time mesh

%output:
%      - g_n : approximated solution of the convolution equation k*g = phi at the point of the mesh t 

%parameters of the mesh
Delta_j = t(2:end)-t(1:end-1);
Delta_max = max(Delta_j);
Nt = length(t)-1;

%parameters for the complex quadrature: see Section 5
Nq = floor(Nt*log(Nt));
[s_l,w_l] = complex_quadrature(1,Nq,Delta_j);

%Algorithm 4

%initialization
g_n = zeros(Nt+1,1);
phi_n = phi(t(2:end));
u_n = zeros(1,Nq);

phi_n = [0 phi_n];
Delta_j = [Delta_max Delta_j];
g_n(2) = phi_n(2)/K_op(1/Delta_j(2));

for n = 3 : Nt+1
    
    %STEP 1
    u_n = u_n./(1-Delta_j(n-1)*s_l) + g_n(n-1)*Delta_j(n-1)./(1-Delta_j(n-1)*s_l);

    %STEP 2
    r_n = phi_n(n) - sum(w_l.*K_op(s_l)./(1-Delta_j(n)*s_l).*u_n);

    %STEP 3
    g_n(n) = r_n/K_op(1/Delta_j(n));    
end

g_n = real(g_n(2:end));