function[Vs] = shotp(f,R,C,I)

% Shot calculat the the shote noise for a wier dipole antenna.
% [Vs] = shot(f,ne,Te,R,C) where;
% f is a vector white the freqensy. 
% ne is the electron densety.
% Te is the electron temperatur.
% R is the antena resistance.
% C is the antena capacitance.
%
% The equastion comes from meyer-vernet artical from 2009(quasi-thermal noise in space plasma) 
%
% shot(f,ne,Te,R,C) deliver a vector white the calulated results of f.
%
% exampel:
%         shot([1 2 3 4 5 6],5e6,1.5e5,5e6,20e-12)
%
% see also nois, ion, electron.
%
% $Id: shot.m,v 1.1 2013/02/26 09:53:00 Pansar Exp $

%% Declaration of constatns (NIST) 


qe= 1.602176565*10^-19; %% Elemetary charge (C)

%% shot noise
for i = f
z= 1./(1./R+(1i.*C.*(2.*pi.*f))); 

Vs= 2.*qe.*I*abs(z).^2./6^2;

end