function phi_n = forward_gcCQ_BDF2(g,K_op,t)

% function that compute the forward generalised convolution quadrature based on the BDF2 
% (for details see [Banjai, Ferrari, Generalized convolution quadrature based on the trapezoidal rule])
% to compute the convolution k * g = phi

%input:
%      - g : given function datum 
%      - K_op : the laplace transform of the kernel k
%      - t : time mesh

%output:
%      - phi_n : approximated convolution k*g = phi at the point of the mesh t 

%parameters of the mesh
Delta_j = t(2:end)-t(1:end-1);
Nt = length(t)-1;

%parameters for the complex quadrature: see Section 5.2
Nq = floor(Nt*log(Nt)^2);
R = 1.5; %See Remark 5
[s_l,w_l] = complex_quadrature(R,Nq,Delta_j);

%Algorithm 1

%initialization
Delta_j = [Delta_j(1) Delta_j];
phi_n = zeros(Nt,1);
g_n = g(t(2:end));
g_n = [0 g_n];

u_n = zeros(1,Nq);
u_n_m_1 = zeros(1,Nq);
u_n_m_2 = zeros(1,Nq);

for n = 3 : Nt+2
    
    %STEP 1 (see also Equation (5.2))
    A_n = Delta_j(n-1)*(Delta_j(n-2)+Delta_j(n-1))/(Delta_j(n-2)+2*Delta_j(n-1));
    B_n = (Delta_j(n-2)+Delta_j(n-1))^2/Delta_j(n-2)/(Delta_j(n-2)+2*Delta_j(n-1));
    C_n = Delta_j(n-1)^2/Delta_j(n-2)/(Delta_j(n-2)+2*Delta_j(n-1));

    u_n = u_n*B_n./(1-A_n*s_l) - u_n_m_2*C_n./(1-A_n*s_l) + g_n(n-1)*A_n./(1-A_n*s_l);
    u_n_m_2 = u_n_m_1;
    u_n_m_1 = u_n;

    %STEP 2
    phi_n(n-2) = sum(w_l.*K_op(s_l).*u_n);

end

phi_n = real(phi_n);