function[Ve] = electron(f,ne,Te,L)

% Electorn calculat the the electron thermal noise for a wier dipole antenna.
% [Ve] = electron(f,ne,Te) where:
% input:
% f is a vector white the freqensy.
% ne is the electron densety.
% Te is the electron temperatur.
% L  antenna lenght
%
% Output:
% Ve is a vector for the voltage power spectrum from the electron nois
%
% The equastion comes from meyer.bernet and chateau artical from 1991(electrostatic noise in non-maxwellian plasmas) 
% 
% electron(f,ne,Te) deliver a vector white the calulated results of f.
%
% exampel:
%         electron([1 2 3 4 5 6],5e6,1.5e5)
%
% see also nois, ion, shot.
%
% $Id: electon.m,v 1.2 2013/04/10 14:23:00 Pansar Exp $

%% Declaration of constatns 

Me= 9.10938291*10^-31;  %% Electron mass                (KG)
KB= 1.3806488*10^-23;   %% Boltsman konstatn            (J/K)
e0= 8.854187817*10^-12; %% Electric constatn            (F/m)
qe= 1.602176565*10^-19; %% Elemetary charge             (C)

k=2;                    %% Kappa

%% Small functions 
v=sqrt((2.*k-3)./k.*KB.*Te./Me);    %% Electron thermisk velocity   (m/s)

Fp=sqrt(ne*qe.^2/(Me*e0))/(2*pi);   %% Plasma frequensy

Ld=sqrt(e0*KB*Te/(ne*(qe)^2));  %% Debay length

n=1;
u=L./Ld;
for i = f 
       
    r= i./Fp;
    
    x=@(z)r.*u./(z.*sqrt(2.*k-1));
    
    F=@(z) 1./x(z).*(sinint(x(z))-1./2.*sinint(2.*x(z))-2./(x(z)).*sin(x(z)./2).^4);
    
    
    el2=@(z)0;
    for j = 0:k
       
        el3=@(z)(factorial(k+j)./factorial(j).*1./((2i).^(k+1+j).*(z+1i).^(k+1-j)));
        
        el2= @(z) el2(z)+el3(z);
    end
    
    el= @(z) 1+z.^2./r.^2.*( 2.*k-1+(-2).^(k+1)./((factorial(factorial(2.*k-3)))).*z.*1i.* el2(z));

  
   int=@(z) z.*(F(z)).*(((1+z.^2).^(k).*abs(el(z)).^2).^(-1));
    
    
   Ve2= integral(int,0,inf);
   
   Ve1=2.^(k+3)./(pi.^2.*e0).*factorial(k)./(factorial(factorial(2.*k-3)).*sqrt(k)).*Me.*v./r.^2;
    
    Ve(n)=Ve1.*Ve2/L^2;   

    n=n+1;
  
end