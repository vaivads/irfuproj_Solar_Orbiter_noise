function[Ve,Vi,Vp,Vt,VE,Va] = noisR(f,Te,Tp,Re,Ri,Rs,I,L,Vb)

% NoisR calculate the thermal noise of ion, electrons and shot 
% for a wire dipole antenna
% [Ve,Vi,Vp,Vs,Vt,VE] = noisR(f,Te,Tp,Re,Ri,Rs,I,L,Vb)
%
% input:
% f is the frequency vector.
% Te is the electron temperatur.
% Tp is the proton temperatur.
% Re plasma resistance of the electron.
% Ri plasma resistance of the ion.
% Rs plasma resistance of the photoelectron.
% I  is the antenna current.
% L  is the antenna lenght.
% Vb volts outside plasma sheet.
%
% Output:
% Ve is a vector for the voltage power spectrum from the electon noise.
% Vi is a vector for the voltage power spectrum from the Ions noise.
% Vp is a vector for the voltage power spectrum from the Photoelctrons
% noise.
% Vs is a vector for the voltage power spectrum from the shote noise.
% VE is a vector for the voltage power spectrum from the e-feild outside
% the plasma sheat noise.
%  
% see also: Noisv11
%
% $Id: nose.m,v 1.1 2013/04/10 14:35:00 Pansar Exp $
%

%%
KB= 1.3800*10^-23;          %% Boltzmann constant           (J/K)

Cp=60e-12;                  %% Plasma capacitance           (F)
Rp=Vb.*L.^2./(4.*KB.*Tp);   %% Impedance of free space      (ohm)

Csc=67e-12;                 %% Spacecraft capacitance       (F)
Rsc=10e3;                   %% Spacecraft resistance        (ohm)

Rin=1e12;                   %% Input resistance             (ohm)
Cin=0;                      %% Input capacitance            (F)

Ra=10e8;                    %% Amplifer resistance          (ohm)
Ca=50e-12;                  %% Amplifer capacitance         (F)

%% tehvene 
TE=sqrt(4*KB*Te*Re);
TI=sqrt(4*KB*Tp*Ri);
TP=sqrt(4*KB*(1/8.61734e-5)*Rs);

n=1;
for w = f
%% Impedance
Zsc=1/(1/Rsc+(Csc*w*2*pi*1i));
Zi=1/(1/Rin+(Cin*w*2*pi*1i));
Za=1/(1/Ra+(Ca*w*2*pi*1i));
ZS=1/(1/Zi);

Zt=Zsc+ZS;
Zp=1/(Cp*w*2*pi*1i)+Rp(n);

%% Electon (Thermal)
SDe=1/(1/Ri+1/Rs+1/Zp+1/Zt)/(Re+1/(1/Ri+1/Rs+1/Zp+1/Zt));

V=TE*abs(SDe);

Ve(n)=V^2/L^2;

%% Ion (Thermal)
SDi=1/(1/Re+1/Rs+1/Zp+1/Zt)/(Ri+1/(1/Re+1/Rs+1/Zp+1/Zt));

V=TI*abs(SDi);

Vi(n)=V^2/L^2;

%% Photoelectron (Thermal)
SDp=1/(1/Re+1/Ri+1/Zp+1/Zt)/(Rs+1/(1/Re+1/Ri+1/Zp+1/Zt));

V=TP*abs(SDp);

Vp(n)=V^2/L^2;


%% Total thermal
if isnan(Vp(n)) == 1
   Vp(n)=inf;
end

Vt(n)=sqrt(Ve(n).^2+Vi(n).^2+Vp(n).^2);

%% QTN
%%SDE=(1/(1/Re+1/Ri+1/Rs+1/Zt))/(1/(1/Re+1/Ri+1/Rs+1/Zt)+Zp); %% serie

SDE=(Zt)/(1/(1/Re+1/Ri+1/Rs+1/Zp)+Zt); %%para

V=sqrt(Vb(n)*L^2)*abs(SDE);
VE(n)=V^2/L^2;

%% Amplyfier
% volt
%SDA=1/(Zsc+1/(1/Re+1/Ri+1/Rs+1/Zp)+1/Za)/(1/(Zsc+1/(1/Re+1/Ri+1/Rs+1/Zp)+1/Za)+Zi);

V=4e-9;
Vav(n)=V^2/L^2;

% current
Z=1/(1/(Zsc+1/(1/Re+1/Ri+1/Rs+1/Zp))+1/(Za)+1/(Zi));
V=10e-15*abs(Z);
Vac(n)=V^2/L^2;

Va(n)=sqrt(Vav(n).^2+Vac(n).^2);
n=n+1;
end

end