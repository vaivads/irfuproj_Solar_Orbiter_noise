function[Vs] = shote(f,ne,T,R,C,L,M)

% Shote calculates the shot noise for a wire dipole antenna.
% [Vs] = shot(f,ne,Te,R,C) where;
% input:
% f is a vector white the frequency. 
% ne is the electron density.
% Te is the electron temperature.
% R is the antenna resistance.
% C is the antenna capacitance.
% L  is the antenna length
% M  is mass
% 
% Output:
% Vs is a vector for the voltage power spectrum from the shot noise
%
% The question comes from Meyer-Vernet article from 2009 (quasi-thermal noise in space plasma) 
%
%
% exampel:
%         shot([1 2 3 4 5 6],5e6,1.5e5,5e6,20e-12)
%
% see also nois, ion, electron.
%
% $Id: shot.m,v 1.2 2013/04/10 14:18:00 Pansar Exp $

%% Declaration of constants (NIST) 

KB= 1.3806488*10^-23;   %% Boltsman constant    J/K)
qe= 1.602176565*10^-19; %% Elemetary charge     (C)

a=0.0116/2;             %% Equivalent radius    (m)
                 
%% small functions

v=sqrt((KB*T)/M);       %%  thermal velocity   (m/s)

Ne=1./sqrt(4.*pi).*ne.*v.*a.*L.*2.*pi;
%% shot noise
n=1;
for w = f
z= 1./(1./R+(1i.*C.*(2.*pi.*w))); 

Vs(n)= 2.*qe.^2*Ne.*abs(z).^2./L^2;
n=n+1;
end
