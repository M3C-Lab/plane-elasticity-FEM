function [j, pqp] = Mapping_lineqp(qp, n1, n2)
% To map quadrature points from [-1, 1] to line(n1, n2)
% Input:
%   qp: The quadrature points sequence in parametric space [-1, 1].
%   n1, n2: The node of line element in the physical space.
% Output:
%   j: The Jacobian for line mapping.
%   pqp: The quadrature points sequence in the physical space.
%       pqp(*, i): The 'i'th point.
%       pqp(1, i): The x-coordinate.
%       pqp(2, i): The y-coordinate.

pqp = zeros(2, length(qp));

for qua = 1 : length(qp)
          pqp(1, qua) = n1(1) + 0.5 * (qp(qua) + 1) * (n2(1) - n1(1));
          pqp(2, qua) = n1(2) + 0.5 * (qp(qua) + 1) * (n2(2) - n1(2));
end

j = 0.5 * sqrt((n1(1)-n2(1))^2 + (n1(2)-n2(2))^2);

end

