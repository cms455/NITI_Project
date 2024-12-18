% LAMB WAVE DISPERSION CURVE EDITOR
% E = Young's Modulus 
% v = Poisson's Ratio
% D = Density
% d = Material Thickness
% mode = Number of modes 
% len = Minimum interval
% maxf = Maximum Frequency.
% vps and vpa are phase velocities of symmetric and antisymmetric modes
% respectively, similarly for vga and vgs
% For Example
%[vps,vpa,vgs,vga] = LambDispersion(7.24e+10,0.33,2780,0.001,8,10,20000)

lambda = 2.37e9;
mu = 26e3;
D = 1000;
d = 710e-6;
mode = 1;
len = 20;
maxf = 7000;

young_modulus = mu*((3*lambda + 2*mu)/(lambda + mu));
poisson_ratio = 0.495;

[vps,vpa,vgs,vga] = LambDispersion(young_modulus,poisson_ratio,D,d,mode,len,maxf);

