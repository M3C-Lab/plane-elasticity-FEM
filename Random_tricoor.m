function point = Random_tricoor()
% To get a random sampling point in a triangle with area coordinates .
% Input:
%   None
% Output:
%   point: The coordinates of the points
%       point(1), points(2), points(3): r, s, t.

p = 100 * rand(3, 1);
sum = p(1) + p(2) + p(3);
point = zeros(3, 1);
point(1) = p(1) / sum;
point(2) = p(2) / sum;
point(3) = p(3) / sum;

end

