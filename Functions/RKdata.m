function [A,b,c,p,q] = RKdata(RK)
%function [A,b,c] = RKdata(RK)
%returns coefficients A,b,c in Butcher notation of Runge-Kutta methods and
% p,q classical order and stage order of the method

%RK = 2x for x-stage Radau IIA (p=2x-1, q=x)
% see [Butcher, Numerical methods for ordinary differential equations, pp. 225-226 (2008)]

%RK = 3x for x-stage Lobatto IIIC pag 228 Butcher2008 (p=2x-2, q=x-1)
% see [Butcher, Numerical methods for ordinary differential equations, pp. 226-229 (2008)]

if (RK==32)

    p=2;q=1;
    A = [1/2 -1/2; 1/2 1/2];
    c = [0 ; 1];
    b = [1/2 ; 1/2];

elseif (RK==33)

    p=4;q=2;
    A = [1/6 -1/3 1/6; 1/6 5/12 -1/12; 1/6 2/3 1/6];
    c = [0 ; 1/2 ; 1];
    b = [1/6 ; 2/3 ; 1/6];

elseif(RK==34)

    p=6;q=3;
    a1 = sqrt(5);
    A = [1/12 -a1/12 a1/12 -1/12; 1/12 1/4 (10-7*a1)/60 a1/60;
        1/12 (10+7*a1)/60 1/4  -a1/60; 1/12 5/12 5/12 1/12];
    b = [1/12; 5/12; 5/12; 1/12];
    c = [0 ;(5-a1)/10; (5+a1)/10; 1];

elseif RK == 35

    p=8;q=4;
    A = [1/20 -7/60 2/15 -7/60 1/20;...
        1/20 29/180 (47-15*sqrt(21))/315 (203-30*sqrt(21))/1260 -3/140;...
        1/20 (329+105*sqrt(21))/2880 73/360 (329-105*sqrt(21))/2880 3/160;...
        1/20 (203+30*sqrt(21))/1260 (47+15*sqrt(21))/315 29/180 -3/140;...
        1/20 49/180 16/45 49/180 1/20];
    c = [0 ; (7-sqrt(21))/14 ; 1/2 ; (7+sqrt(21))/14; 1];
    b = [1/20 ; 49/180 ; 16/45 ; 49/180 ; 1/20];

elseif (RK==22)

    p=3;q=2;
    A = [5/12 -1/12; 3/4 1/4]; c = [1/3; 1];
    b = [3/4; 1/4];

elseif (RK==23)

    p=5;q=3;
    A = [(88-7*sqrt(6))/360 (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225;
        (296+169*sqrt(6))/1800 (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225;
        (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
    b = [(16-sqrt(6))/36 (16+sqrt(6))/36 1/9].';
    c = [(4-sqrt(6))/10 ; (4+sqrt(6))/10 ; 1];

end