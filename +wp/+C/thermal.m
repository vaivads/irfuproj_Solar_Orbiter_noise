function[Vt] = thermal(f,Re,Rp,Rpe,Rt,C,Te,Tp)

% Shot calculat the the shote noise for a wier dipole antenna.
% [Vs] = shot(f,ne,Te,R,C) where;
% f is a vector white the freqensy. 
% ne is the electron densety.
% Te is the electron temperatur.
% R is the antena resistance.
% C is the antena capacitance.
%
% The equastion comes from ????? 
%
% shot(f,ne,Te,R,C) deliver a vector white the calulated results of f.
%
% exampel:
%         shot([1 2 3 4 5 6],5e6,1.5e5,5e6,20e-12)
%
% see also nois, ion, electron.
%
% $Id: shot.m,v 1.1 2013/02/26 09:53:00 Pansar Exp $


%% Declaration of constatns 
qe= 1.6022*10^-19;      %% Elemetary charge             (C)
K= 8.61734e-5;          %% Boltsman konstatn for temp   (eV/K)


%% shot noise
for i = f
R1=Re; 
R2=Rp;
R3=Rpe;
%Rp=R1.*R2./(R1+R2);
%Rt=Rp.*R3./(Rp+R3);
Z=1./(1./Rt+(1i.*C.*(2.*pi.*f)));

Vt= 4.*qe.*abs(Z).^2.*((K.*Te)./R1+(K.*Tp)./R2+1./R3)./(6^2);

end
