clc;clear all;
addpath('Preprocess', 'Mapping & Quadrature', 'Postprocess')

[qp,w] = Gauss(5,-1,1);
n1 = [2,2];n2 = [3,4];
[j,pqp] = Mapping_lineqp(qp, n1, n2);