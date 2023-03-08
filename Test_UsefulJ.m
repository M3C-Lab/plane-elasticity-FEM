clc; clear all;

% Physical space
n1 = [0, 0];
n2 = [2, 1];
n3 = [0, 2];

Jacobian_matrix = [1, n1(1), n1(2);
                   1, n2(1), n2(2);
                   1, n3(1), n3(2)];
J = det(Jacobian_matrix);

Quad_degree = 4;

[qp, wq, nqp] = TriangularQuad(Quad_degree);

f1 = @(x, y) x * y;
% ^f1(r,s,t) = f1(n1(1)*r + n2(1)*s + n3(1)*t, n1(2)*r + n2(2)*s + n3(2)*t
f2 = @(x, y) x^2 + y^2;
f3 = @(x, y) x^3 + y^3;

result1 = 0;
result2 = 0;
result3 = 0;

for qua = 1 : nqp
    result1 = result1 + J * wq(qua) * ...
        f1(n1(1)*qp(1, qua) + n2(1)*qp(2, qua) + n3(1)*qp(3, qua),...
           n1(2)*qp(1, qua) + n2(2)*qp(2, qua) + n3(2)*qp(3, qua));
    result2 = result2 + J * wq(qua) * ...
        f2(n1(1)*qp(1, qua) + n2(1)*qp(2, qua) + n3(1)*qp(3, qua),...
           n1(2)*qp(1, qua) + n2(2)*qp(2, qua) + n3(2)*qp(3, qua));
    result3 = result3 + J * wq(qua) * ...
        f3(n1(1)*qp(1, qua) + n2(1)*qp(2, qua) + n3(1)*qp(3, qua),...
           n1(2)*qp(1, qua) + n2(2)*qp(2, qua) + n3(2)*qp(3, qua));
end
