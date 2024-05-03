function g_n = backward_gcCQ_BDF2(phi,K_op,t)

% function that compute the backward generalised convolution quadrature based on the BDF2 
% (for details see [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule])
% to solve convolution equation k * g = phi

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

%parameters for the complex quadrature: see Section 5.2
Nq = floor(Nt*log(Nt)^2);
R = 1.5; %See Remark 5
[s_l,w_l] = complex_quadrature(R,Nq,Delta_j);

%Algorithm 2

%initialization
g_n = zeros(Nt+2,1);
phi_n = phi(t(2:end));
u_n = zeros(1,Nq);
u_n_m_1 = zeros(1,Nq);
u_n_m_2 = zeros(1,Nq);

phi_n = [0 0 phi_n];
Delta_j = [Delta_max Delta_max Delta_j];
g_n(3) = phi_n(3)/K_op(3/2/Delta_j(3));

for n = 4 : Nt+2
    
    %STEP 1 (see also Equation (5.2))
    A_n_m_1 = Delta_j(n-1)*(Delta_j(n-2)+Delta_j(n-1))/(Delta_j(n-2)+2*Delta_j(n-1));
    B_n_m_1 = (Delta_j(n-2)+Delta_j(n-1))^2/Delta_j(n-2)/(Delta_j(n-2)+2*Delta_j(n-1));
    C_n_m_1 = Delta_j(n-1)^2/Delta_j(n-2)/(Delta_j(n-2)+2*Delta_j(n-1));
    
    u_n = u_n*B_n_m_1./(1-A_n_m_1*s_l) - u_n_m_2*C_n_m_1./(1-A_n_m_1*s_l) + g_n(n-1)*A_n_m_1./(1-A_n_m_1*s_l);
    u_n_m_2 = u_n_m_1;
    u_n_m_1 = u_n;

    %STEP 2
    A_n = Delta_j(n)*(Delta_j(n-1)+Delta_j(n))/(Delta_j(n-1)+2*Delta_j(n));
    B_n = (Delta_j(n-1)+Delta_j(n))^2./Delta_j(n-1)/(Delta_j(n-1)+2*Delta_j(n));
    C_n = Delta_j(n)^2/Delta_j(n-1)/(Delta_j(n-1)+2*Delta_j(n));
    
    r_n = phi_n(n) - sum(w_l.*K_op(s_l).*(u_n.*B_n./(1-A_n*s_l) - u_n_m_2.*C_n./(1-A_n*s_l)));
    
    %STEP 3
    g_n(n) = r_n/K_op(1/A_n);
    
end

g_n = real(g_n(3:end));


