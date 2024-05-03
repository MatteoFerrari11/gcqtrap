function g_n = backward_gcCQ_RK(phi,K_op,t,RK)

% function that compute the backward generalised convolution quadrature
% based a stiffly accurate Runge-Kutta method
% (for details see [Lopez Fernandez, Sauter, Generalized convolution quadrature based on Runge-Kutta methods].
% and [Leitner, Schanz, Generalized convolution quadrature based boundary element method for uncoupled 
% thermoelasticity (Section 3)])

%input:
%      - phi : given function datum 
%      - K_op : the laplace transform of the kernel k
%      - t : time mesh
%      - RK : the Runge-Kutta method

        %RK = 2x for x-stage Radau IIA
        %RK = 3x for x-stage Lobatto IIIC

%output:
%      - g : approximated solution of the convolution equation k*g = phi at the point of the mesh t 

[A,b,c] = RKdata(RK);
s = length(b); %number of stages

es = [zeros(s-1,1);1];

%parameters of the mesh
Delta_j = t(2:end)-t(1:end-1);
Nt = length(t)-1;

%parameters for the complex quadrature: see Section 5.2
Nq = 2*floor(Nt*log(Nt)^2);
R = 5;
[s_l,w_l] = complex_quadrature(R,Nq,Delta_j);

%initialization
g_n = zeros(s,Nt);
phi_n = zeros(s,Nt);
for i = 1 : s
    for n = 1 : Nt
        phi_n(i,n) = phi(t(n)+c(i)*Delta_j(n));
    end
end
u_n = zeros(s,Nq);

[V,D] = eig(A*Delta_j(1));
ev = K_op(1./diag(D));
K1 = V*diag(ev)/V;
g_n(:,1) = K1\phi_n(:,1);

for n = 2 : Nt
    
    %STEP 1
    for l = 1 : Nq
        u_n(:,l) = (eye(size(A))-Delta_j(n-1)*s_l(l)*A)\((kron(ones(s,1),es.')*u_n(:,l)+Delta_j(n-1)*A*g_n(:,n-1)));
    end
    
    %STEP 2
    r_n = phi_n(:,n);
    for l = 1 : Nq
        r_n = r_n - w_l(l)*K_op(s_l(l))*u_n(end,l)*inv(eye(size(A))-Delta_j(n)*s_l(l)*A)*ones(s,1);
    end
    
    %STEP 3
    [V,D] = eig(A*Delta_j(n));
    ev = K_op(1./diag(D));
    Kn = V*diag(ev)/V;
    g_n(:,n) = Kn\r_n;
    
end

g_n = real(g_n(end,:)).';