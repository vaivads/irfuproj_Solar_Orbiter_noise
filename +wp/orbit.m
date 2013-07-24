function[] = orbit()
%Orbit calculate the noise magnitude over almost the full travel distance
%of the Solar Orbiter and is limited to a case where the antenna is 
%resistive coupled to the plasma.
%
%[] = orbit() do not accepts any input valuse. Changes must be done within
%the program.
%
% The result will be plotted in a semilogarithmic scale.
%
% IT is the vector for the biased current over the distance to sun
%
% see also data2, shot, ion, electron, noisR.
%
% $Id: nose.m,v 1.1 2013/04/10 14:52:00 Pansar Exp $

load 'wp/Helios_data'

Units=irf_units;
Me= Units.me;      %% Electron mass           (Kg)
L=5;               %% Antenna lenght          (m)
r=0.575e-2;        %% Antenna radius          (m)

D=0.30:0.02:0.9;
P=[50 95 95];
Vsweep=[-5:0.00001:15];
IT=[-16e-6 -13e-6 -12e-6 -10e-6 -9e-6 -8e-6 -7e-6 -6.5e-6 -6e-6 -5.5e-6...
    -5e-6 -4.5e-6 -4e-6 -3.5e-6 -3e-6 -3e-6 -3e-6 -2.5e-6 -2.5e-6 -2e-6...
    -2e-6 -2e-6 -1.5e-6 -1.5e-6 -1.5e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6];
f=10.^[log10(10^0):2:log10(1e1)];

for j = 1:length(D)
[Ne,Tp,B,V,RBt,RBe,RBi,RBp,RUt,RUe,RUi,RUp,VU,VB,Is,Ib,Iu,IT0] = wp.data2(D(j),IT(j),P,Vsweep,r);

Ne=Ne.*1e6;         %% Electron density             (m^-3) 
Te=Tp./2.5;         %% Electron temperatur          (K)
V=V.*1e3;           %% Solar wind velocity          (m/s)
B=B.*1e-9;          %% Magnetic field               (T)

Cp=30e-12;          %% Plasma capacitance           (F)

%% Calculations
for i=1:length(P),

I{i}=wp.C.ion(f,Ne(i),Te(i),Tp(i),V(i),L); % QTN ion

SBe{i}=wp.C.shote(f,Ne(i),Te(i),RBt(i),Cp,L,Me); %shot noise

SUe{i}=wp.C.shote(f,Ne(i),Te(i),RUt(i),Cp,L,Me); % shot noise

%% Biased 
IB(i)=sqrt(Is(i)^2+(Is(i)-Ib(i))^2); 
IU(i)=sqrt(Is(i)^2+(Is(i)-Iu(i))^2);

% Other noise
[Veb{i},Vib{i},Vpb{i},Vb{i},VEb{i},Vab{i}] = wp.R.noisR(f,Te(i),Tp(i),RBe(i),RBi(i),RBp(i),IB(i),L,I{i});

[Veu{i},Viu{i},Vpu{i},Vu{i},VEu{i},Vau{i}] = wp.R.noisR(f,Te(i),Tp(i),RUe(i),RUi(i),RUp(i),IU(i),L,I{i});

Vb{i}=Vb{i}+VEu{i}+Vau{i}+SBe{i};
Vu{i}=Vu{i}+VEu{i}+Vau{i}+SUe{i};
end
% Calulation a mean value
T1(j)=mean(Vb{1});
T2(j)=mean(Vb{2});
T3(j)=mean(Vb{3});
U1(j)=mean(Vu{1});
U2(j)=mean(Vu{2});
U3(j)=mean(Vu{3});

NE1(j)=mean(Ne(1));
NE2(j)=mean(Ne(2));
NE3(j)=mean(Ne(3));

TE1(j)=mean(Te(1));
TE2(j)=mean(Te(2));
TE3(j)=mean(Te(3));

VE1(j)=mean(V(1));
VE2(j)=mean(V(2));
VE3(j)=mean(V(3));

BE1(j)=mean(B(1));
BE2(j)=mean(B(2));
BE3(j)=mean(B(3));

Ig(j)=mean(Is);
IT1(j)=mean(IT0(1));
IT2(j)=mean(IT0(2));
IT3(j)=mean(IT0(3));

j/length(D)
 
end

% Plot over the distance for a biased antenna           
Bias=figure(1);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

semilogy(D,T1,'-b',D,T2,'-k',D,T3,'-r')
legend('Average','High particle temperature','High particle density ','Location','Best')

xlabel('Distance to the sun (AU)')
ylabel('Total noise((V/m)^2/hz)')

set(Bias,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Biasdis','.eps']; %location of saved
% plot
print( '-depsc2' , name ) % saving plot

% plot of the proton temperature over the distance
Temp=figure(2);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

semilogy(D,TE1,'-b',D,TE2,'-k',D,TE3,'-r')
legend('Average','High particle temperature','High particle density ','Location','Best')
xlabel('Distance to the sun (AU)')
ylabel('Particle temperature (K)')
set(Temp,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Temp','.eps']; %location of saved
%plot
print( '-depsc2' , name ) %saving plot

% plot of the particle density over the distance
Dens=figure(3);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

semilogy(D,NE1,'-b',D,NE2,'-k',D,NE3,'-r')
legend('Average','High particle temperature','High particle density ','Location','Best')
xlabel('Distance to the sun (AU)')
ylabel('Particle density (m^3)')
set(Dens,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Dens','.eps'];
print( '-depsc2' , name )

% plot of the solar wind speed over the distance
Wind=figure(4);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

plot(D,VE1,'-b',D,VE2,'-k',D,VE3,'-r')
legend('Average','High particle temperature','High particle density ','Location','Best')
xlabel('Distance to the sun (AU)')
ylabel('Solar wind speed (m/s)')
set(Wind,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Wind','.eps']; %location of saved
%plot
print( '-depsc2' , name ) %saving plot

% plot of the magnetic field over the distance
Mag=figure(5);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

semilogy(D,BE1,'-b',D,BE2,'-k',D,BE3,'-r')
legend('Average','High particle temperature','High particle density ','Location','Best')
xlabel('Distance to the sun (AU)')
ylabel('Magentic field (T)')
set(Mag,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Mag','.eps'];%location of saved
%plot
print( '-depsc2' , name )% saving plot

% plot of the current inside the antenna
Current=figure(6);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

plot(D,Ig,'-m',D,IT,'-c',D,IT1,'-b',D,IT2,'-k',D,IT3,'-r')
xlabel('Distance to the sun (AU)')
ylabel('Current (A)')
legend('Photoelectron saturation','Biased','Average','High particle temperature','High particle density ','Location','Best')
set(Current,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Current','.eps'];%location of saved
%plot
print( '-depsc2' , name ) %saving plot

% plot over the distance for a unbiased antenna
Unbias=figure(7);
% setup of plot
set(0,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'PaperUnits','centimeters')
xSize = 24; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

semilogy(D,U1,'-b',D,U2,'-k',D,U3,'-r')
legend('Average','High particle temperature','High particle density ','Location','Best')
xlabel('Distance to the sun (AU)')
ylabel('Total noise ((V/m)^2/Hz)')
set(Unbias,'color','white'); % white background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Unbiasdis','.eps'];%location of saved
%plot
print( '-depsc2' , name ) %saving plot
end