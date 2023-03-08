clc; clear all;

% Parent space
p1 = [0, 0];
p2 = [1, 0];
p3 = [0, 1];

% Physical space
n1 = [2, 2];
n2 = [0, 1.5];
n3 = [0.5, 0];

[para2phys,phys2para] = Mapping(n1, n2, n3, p1, p2, p3);
%       x = para2phys(1,1) * xi + para2phys(1,2) * eta + para2phys(1,3)
%       y = para2phys(2,1) * xi + para2phys(2,2) * eta + para2phys(2,3)

Jacobian_matrix = [para2phys(1, 1), para2phys(1, 2);
                   para2phys(2, 1), para2phys(2, 2)];
J = det(Jacobian_matrix);

Quad_degree = 5;

[qp_p, wq_p, nqp_p] = TriangularQuad(Quad_degree, p1, p2, p3);
[qp_n, wq_n, nqp_n] = TriangularQuad(Quad_degree, n1, n2, n3);

f1 = @(x, y) x * y;
% ^f1(xi, eta) = f1(para2phys(1, 1)*xi + para2phys(1, 2)*eta+ para2phys(1, 3),
%                  para2phys(2, 1)*xi + para2phys(2, 2)*eta+ para2phys(2, 3))
f2 = @(x, y) x^2 + 50 * y^2;
f3 = @(x, y) x^3 + 100 * y^4;

result1_p = 0;
result1_n = 0;
result2_p = 0;
result2_n = 0;
result3_p = 0;
result3_n = 0;


for qua = 1 : nqp_p
    result1_p = result1_p + wq_p(qua) * ...
        f1(para2phys(1, 1)*qp_p(1,qua) + para2phys(1, 2)*qp_p(2,qua)+ para2phys(1, 3),...
           para2phys(2, 1)*qp_p(1,qua) + para2phys(2, 2)*qp_p(2,qua)+ para2phys(2, 3));
    result2_p = result2_p + wq_p(qua) * ...
        f2(para2phys(1, 1)*qp_p(1,qua) + para2phys(1, 2)*qp_p(2,qua)+ para2phys(1, 3),...
           para2phys(2, 1)*qp_p(1,qua) + para2phys(2, 2)*qp_p(2,qua)+ para2phys(2, 3));
    result3_p = result3_p + wq_p(qua) * ...
        f3(para2phys(1, 1)*qp_p(1,qua) + para2phys(1, 2)*qp_p(2,qua)+ para2phys(1, 3),...
           para2phys(2, 1)*qp_p(1,qua) + para2phys(2, 2)*qp_p(2,qua)+ para2phys(2, 3)); 
end

for qua = 1 : nqp_n
    result1_n = result1_n + wq_n(qua) * f1(qp_n(1,qua), qp_n(2,qua));
    result2_n = result2_n + wq_n(qua) * f2(qp_n(1,qua), qp_n(2,qua));
    result3_n = result3_n + wq_n(qua) * f3(qp_n(1,qua), qp_n(2,qua));
end

