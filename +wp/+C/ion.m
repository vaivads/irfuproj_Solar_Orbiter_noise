function[Vi] = ion(f,ne,Te,Tp,V,L)

% Ion calculat the thermal ion thermal noise for a wier dipole antenna.
% [Vi] = ion(f,ne,Te,Tp,V) where:
% input:
% F is a vector white the freqensy.
% ne is the electron densety.
% Te is the electron temperatur.
% Tp is the proten temperatur.
% V is the solar wind speed.
%
% Output:
% Vi is a vector for the voltage power spectrum from the ion nois
%
% The equastion comes from issautier artical from 1999(quasi-thermal noise in drifting plasma) 
%
% ion(f,ne,Te,Tp,V) deliver a vector white the calulated results of f.
%
% exampel:
%         ion([1 2 3 4 5 6],5e6,1.5e5,0.8e5,354e3)
%
% see also nois, shot, electron.
%
% $Id: ion.m,v 1.2 2013/04/10 14:21:00 Pansar Exp $


%% Declaration of constatns

Me=9.1094e-31;      %% Electron mass                (KG)
KB=1.3800e-23;      %% Boltsman konstatn            (J/K)
e0=8.8542e-12;      %% Electric constatn            (F/m)
qe=1.6022e-19;      %% Elemetary charge             (C)

%% Small functions 
n=1;
v=sqrt((KB*Te)/Me);             %% Electron thermisk velocity   (m/s)
Ld=sqrt(e0*KB*Te/(ne*(qe)^2));  %% Debay length

t= Te/Tp;
M= V/v;
u= L/Ld;

for i = f 
    
OM=2.*pi.*i.*Ld./V;

%% Proton thermal noise  

IPt= @(y,s) y.*(sin(s./2).^4./(s.^2.*((y.*u).^2-s.^2).^(1./2)))./((y.^2+1+OM.^2).*(y.^2+1+OM.^2+t));

V2=integral2(IPt,0,inf,0,@(y) y.*u);

V1= 64/pi*sqrt(2*Me*KB)/(4*pi*e0)*sqrt(Te)/M;

Vi(n)=(V1.*V2)/L^2;
n=n+1;
end

end