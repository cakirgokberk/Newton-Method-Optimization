% Optimization Theory
% Class Assignment -  Newton Algorithm Main File
% 28/11/2020 

clc;
clear;
close all;

f = @(x1,x2) (100*(x2 - x1.^2).^2 + (1-x1).^2);

gradient = @(x1,x2)  [(100*(-4*x2*x1 + 4*x1^3)+(-2 +2*x1));(100*(2*x2-2*x1^2))] ;
hessian  = @(x1,x2)  [100*(-4*x2 +12*x1^2)+2,-400*x1;-400*x1,200];
x0 = [-2 -2]';

lowerLimit = -3;
upperLimit =  3;

e1 = 1e-9;
e2 = 1e-9;
e3 = 1e-9;

N_Max = 200;
tic
[Y,X,k] = alteredNewton(f,hessian,gradient,N_Max,x0,e1,e2,e3,lowerLimit,upperLimit)
toc

figure (1)
ezsurf(f);



