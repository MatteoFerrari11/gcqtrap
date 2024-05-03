function phi_n = forward_gcCQ_RK(g,K_op,t,RK)

% function that compute the forward generalised convolution quadrature
% based a stiffly accurate Runge-Kutta method
% (for details see [Lopez Fernandez, Sauter, Generalized convolution quadrature based on Runge-Kutta methods].
% and [Leitner, Schanz, Generalized convolution quadrature based boundary element method for uncoupled
% thermoelasticity (Section 3)]) to compute the convolution k * g = phi

%input:
%      - g : given function datum
%      - K_op : the laplace transform of the kernel k
%      - t : time mesh
%      - RK : the Runge-Kutta method

%RK = 2x for x-stage Radau IIA
%RK = 3x for x-stage Lobatto IIIC

%output:
%      - phi_n : approximated convolution k*g = phi at the point of the mesh t_n

[A,b,c] = RKdata(RK);
s = length(b); %number of stages

%parameters of the mesh
Delta_j = t(2:end)-t(1:end-1);
Nt = length(t)-1;

%parameters for the complex quadrature: see Section 5.2
Nq = 2*floor(Nt*log(Nt)^2);
R = 5;
[s_l,w_l] = complex_quadrature(R,Nq,Delta_j);

g_n = zeros(s,Nt);
for i = 1: s
    for n = 1 : Nt
        g_n(i,n) = g(t(n)+c(i)*Delta_j(n));
    end
end

g = zeros(s,Nt);
Un = zeros(s,Nq);
un = zeros(1,Nq);

for n = 1 : Nt

    for l = 1 : Nq
        Un(:,l) = (eye(size(A))-Delta_j(n)*s_l(l)*A)\(un(l)*ones(size(b))+Delta_j(n)*A*g_n(:,n));
    end
    un = un.*(1-b'*inv(A)*ones(size(b))) + b'*inv(A)*Un;

    for l = 1 : Nq
        g(:,n) = g(:,n) + w_l(l)*K_op(s_l(l))*Un(:,l);
    end
end


phi_n(:) = g(end,:);

phi_n = real(phi_n)';