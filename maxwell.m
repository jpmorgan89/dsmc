function [u, v, w] = maxwell(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R1 = rand([n 1]);
R2 = rand([n 1]);

theta = 2*pi*R1;
r = (-log(R2)).^(1/2);

u = r.*cos(theta);
v = r.*sin(theta);

R1 = rand([n 1]);
R2 = rand([n 1]);

theta = 2*pi*R1;
r = (-log(R2)).^(1/2);

w = r.*cos(theta);
end

