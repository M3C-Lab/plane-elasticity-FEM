function [qp, wq, nqp] = TriangularQuad(degree, p1, p2, p3)
% To generate Hammer's quadrature info on Page 173 for the default 
% parent triangle on Page 121 with parent coordinates (xi, eta).
% Input:
%   degree: The degree of precision.
%   p1,p2,p3: The coordinates of the nodes of the parent triangle.
% Output:
%   qp: The coodinates in the parameter space (xi, eta).
%       qp(*, e): The 'e'th quadrature points.
%       qp(1, *): The 'xi' coordinate.
%       qp(2, *): The 'eta' coordinate.
%   wq: The quadrature weight.
%       wq(e): The weight of the 'e'th quadrature point.
%   nqp: The number of quadrature points.

% !!!!!! Note: We will make the number of quadrature points (nqp) as little 
%   as possible. The available choices are:
%   (nqp, degree) = 
%   (3, 2)  (4, 3)  (6, 4)  (7, 5)  (12, 6) (13, 7)
% !!!!!! Note:  The triangular coordinates (r, s, t) would be mapped into
%   parameter coodinates(xi, eta)

if degree == 2
    nqp = 3;
    triangular_qp = zeros(3, nqp); % qp with triangular coordinates
    wq = zeros(1, nqp);
    for ii = 1 : 3
        for jj = 1 : 3
            if ii == jj
                triangular_qp(ii, jj) = 0.666666666666667;
            else
                triangular_qp(ii, jj) = 0.166666666666667;
            end
            wq(1, jj) = 0.333333333333333;
        end
    end
    
elseif degree == 3
    nqp = 4;
    triangular_qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    for ii = 1 : 3
        triangular_qp(ii, 1) = 0.333333333333333;
        for jj = 2 : 4
            if ii == jj - 1
                triangular_qp(ii, jj) = 0.6;
            else
                triangular_qp(ii, jj) = 0.2;
            end
            wq(1, jj) = 0.520833333333333;
        end
    end
    wq(1, 1) = -0.56250;

elseif degree == 4
    nqp = 6;
    triangular_qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    for ii = 1 : 3
        for jj = 1 : 3
            wq(1, jj) = 0.109951743655322;
            if ii == jj
                triangular_qp(ii, jj) = 0.816847572980459;
            else
                triangular_qp(ii, jj) = 0.091576213509771;
            end
        end
        for jj = 4 : 6
            wq(1, jj) = 0.223381589678011;
            if ii == jj - 3
                triangular_qp(ii, jj) = 0.108103018168070;
            else
                triangular_qp(ii, jj) = 0.445948490915965;
            end
        end
    end
    
elseif degree == 5
    nqp = 7;
    triangular_qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    wq(1, 1) = 0.22500;
    for ii = 1 : 3
        triangular_qp(ii, 1) = 0.333333333333333;
        for jj = 2 : 4
            if ii == jj - 1
                triangular_qp(ii, jj) = 0.797426985353087;
            else
                triangular_qp(ii, jj) = 0.101286507323456;
            end
            wq(1, jj) = 0.125939180544827;
        end
        for jj = 5 : 7
            if ii == jj - 4
                triangular_qp(ii, jj) = 0.059715871789770;
            else
                triangular_qp(ii, jj) = 0.470142064105115;
            end
            wq(1, jj) = 0.132394152788506;
        end
    end
    
elseif degree == 6
    nqp = 12;
    triangular_qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    for ii = 1 : 3
        for jj = 1 : 3
            wq(1, jj) = 0.050844906370207;
            if ii == jj
                triangular_qp(ii, jj) = 0.873821971016996;
            else
                triangular_qp(ii, jj) = 0.063089014491502;
            end
        end
        for jj = 4 : 6
            wq(1, jj) = 0.116786275726379;
            if ii == jj - 3
                triangular_qp(ii, jj) = 0.501426509658179;
            else
                triangular_qp(ii, jj) = 0.249286745170910;
            end
        end
        for jj = 7 : 12
            wq(1, jj) = 0.082851075618374;
        end
        r = 0.636502499121399;
        s = 0.310352451033785;
        t = 0.053145049844816;
        triangular_qp(:, 7) = [r, s, t]';
        triangular_qp(:, 8) = [s, t, r]';
        triangular_qp(:, 9) = [t, r, s]';
        triangular_qp(:, 10) = [r, t, s]';
        triangular_qp(:, 11) = [t, s, r]';
        triangular_qp(:, 12) = [s, r, t]';
    end
    
elseif degree == 7
    nqp = 13;
    triangular_qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    wq(1, 1) = -0.149570044467670;
    for ii = 1 : 3
        triangular_qp(ii, 1) = 0.333333333333333;
        for jj = 2 : 4
            wq(1, jj) = 0.175615257433204;
            if ii == jj - 1
                triangular_qp(ii, jj) = 0.479308067841923;
            else
                triangular_qp(ii, jj) = 0.260345966079038;
            end
        end
        for jj = 5 : 7
            wq(1, jj) = 0.053347235608839;
            if ii == jj - 4
                triangular_qp(ii, jj) = 0.869739794195568;
            else
                triangular_qp(ii, jj) = 0.065130102902216;
            end
        end
        for jj = 8 : 13
            wq(1, jj) = 0.077113760890257;
        end
        r = 0.638444188569809;
        s = 0.312865496004875;
        t = 0.086903154253160;
        triangular_qp(:, 7) = [r, s, t]';
        triangular_qp(:, 8) = [s, t, r]';
        triangular_qp(:, 9) = [t, r, s]';
        triangular_qp(:, 10) = [r, t, s]';
        triangular_qp(:, 11) = [t, s, r]';
        triangular_qp(:, 12) = [s, r, t]';
    end
    
else
        disp('TriangularQuad: Please check the data input.');
        disp("The degree input should be a interger from 2 to 7.");
        return;
end

% Map qp into (xi,eta) space.
qp = zeros(2, nqp);

for kk = 1 : nqp
    qp(1, kk) = p1(1) * triangular_qp(1, kk) + p2(1) * triangular_qp(2, kk)...
        + p3(1) * triangular_qp(3, kk);
    qp(2, kk) = p1(2) * triangular_qp(1, kk) + p2(2) * triangular_qp(2, kk)...
        + p3(2) * triangular_qp(3, kk);
end

end

