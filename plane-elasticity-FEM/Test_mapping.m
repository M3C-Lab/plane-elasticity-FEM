clc;clear all;
addpath('Preprocess', 'Mapping & Quadrature', 'Postprocess')

p1 = [2, 2]';
p2 = [0, 1.5]';
p3 = [0.5, 0]';

[J, matrix] = Mapping_p2rst(p1, p2, p3);

inver = inv(matrix);