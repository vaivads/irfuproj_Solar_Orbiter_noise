function[Vs] = shote(f,ne,T,R,C,L,M)

% Shot calculat the thermal shote noise for a wier dipole antenna.
% [Vs] = shot(f,ne,Te,R,C) where;
% input:
% f is a vector white the freqensy. 
% ne is the electron densety.
% Te is the electron temperatur.
% R is the antena resistance.
% C is the antena capacitance.
% 
% Output:
% Vs is a vector for the voltage power spectrum from the shot nois
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
% $Id: shot.m,v 1.2 2013/04/10 14:18:00 Pansar Exp $

%% Declaration of constatns (NIST) 

KB= 1.3806488*10^-23;   %% Boltsman konstatn            (J/K)
qe= 1.602176565*10^-19; %% Elemetary charge             (C)

a=0.0116/2;             %% equivalent radius            (m)
                 
%% small funtions

v=sqrt((KB*T)/M);           %%  thermisk velocity   (m/s)

Ne=1./sqrt(4.*pi).*ne.*v.*a.*L.*2.*pi;
%% shot noise
n=1;
for w = f
z= 1./(1./R+(1i.*C.*(2.*pi.*w))); 

Vs(n)= 2.*qe.^2*Ne.*abs(z).^2./L^2;
n=n+1;
end