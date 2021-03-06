function[Vi] = ion(f,ne,Te,Tp,V,L)

% Ion calculates the ion QTN noise for a wire dipole antenna.
% [Vi] = ion(f,ne,Te,Tp,V) where:
% input:
% F is the frequency vector
% ne is the ion density.
% Te is the electron temperature.
% Tp is the Proton Temperatur.
% V is the solar wind speed.
% L  is the antenna length
%
% Output:
% Vi is a vector for the voltage power spectrum of the ion QTN noise
%
% The equastion is based upon the Issautires paper "Quasi-thermal
% noise in drifting plasma" from 1999.
%
% exampel:
%         ion([1 2 3 4 5 6],5e6,1.5e5,0.8e5,354e3,5)
%
% see also nois, shot, electron.
%
% $Id: ion.m,v 1.2 2013/04/10 14:21:00 Pansar Exp $


%% Declaration of constants

Me=9.1094e-31;      %% Electron mass                (KG)
KB=1.3800e-23;      %% Boltzmann constant           (J/K)
e0=8.8542e-12;      %% Electric constant            (F/m)
qe=1.6022e-19;      %% Elementary charge            (C)

%% Small functions 
n=1;
v=sqrt((KB*Te)/Me);             %% Electron thermal velocity   (m/s)
Ld=sqrt(e0*KB*Te/(ne*(qe)^2));  %% Debye length

t= Te/Tp;
M= V/v;
u= L/Ld;

for i = f 
    
OM=2.*pi.*i.*Ld./V;

%% Proton QTN noise  

IPt= @(y,s) y.*(sin(s./2).^4./(s.^2.*((y.*u).^2-s.^2).^(1./2)))./((y.^2+1+OM.^2).*(y.^2+1+OM.^2+t));

V2=integral2(IPt,0,inf,0,@(y) y.*u);

V1= 64/pi*sqrt(2*Me*KB)/(4*pi*e0)*sqrt(Te)/M;

Vi(n)=(V1.*V2)/L^2;
n=n+1;
end

end
