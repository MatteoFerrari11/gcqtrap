function [s_l,w_l] = complex_quadrature(R,Nq,Delta_j)

% function that compute the parameters for complex quadrature required for gCQ methods
% (for details see
% -[Lopez Fernandez, Sauter, Fast and stable contour integration for high order 
%       divided differences via elliptic functions]
% -[Lopez Fernandez, Sauter, Generalized convolution quadrature with variable time
%   stepping. Part II: Algorithm and numerical results]
% -[Banjai Ferrari, Generalized convolution quadrature based on the
%   trapezoidal rule (Section 5.2)]

%this function requires Driscoll's toolbox on elliptic functions

%input:
%      - R : this parameter modifies the radius of the complex circle
%      used for complex integration. 
%      Suggested choice are:
%           R = 1 for BDF1
%           R = 1.5 for BDF2
%           R = 2 for trapezoidal rule
%           R = 5 for Radau IIA
%      See [Banjai Ferrari, Generalized convolution quadrature based on the
%      trapezoidal rule (Remark 4)]

%      - Nq : number of quadrature points used in the the complex integration.
%      Suggested choice are:
%           Nq = Nt*log(Nt) for BDF1
%           Nq = Nt*log(Nt)^2 for BDF2
%           Nq = Nt*log(Nt)^2 for trapezoidal rule
%           Nq = 2*Nt*log(Nt)^2 for Radau IIA
%       where Nt is the number of time steps in the gCQ      

%      - Delta_j : distance between consecutive time steps in the gCQ

%output:
%      - s_l : nodes   for the complex quadrature
%      - w_l : weights for the complex quadrature

Delta_min = min(Delta_j);
Delta_max = max(Delta_j);

M = R*max(Delta_max^(-2),Delta_min^(-1)); 
q = M*Delta_max;
k = (q-sqrt(2*q-1))/(q+sqrt(2*q-1));

[K,d_K] = ellipkkp(-log(k^2)/pi/2);
sigma_l = -K + ((1:Nq)-1/2)*4*K/Nq + 1i/2*d_K;
[sn,cn,dn] = ellipjc(sigma_l,-log(k^2)/pi/2);

d_gamma = M*sqrt(2*q-1)/(q-1)*2*cn.*dn./k./(k^(-1)-sn).^2;
gamma = M/(q-1)*(sqrt(2*q-1)*(k^(-1)+sn)./(k^(-1)-sn)-1);

s_l = gamma;
w_l = 4*K/2/pi/1i*d_gamma/Nq;