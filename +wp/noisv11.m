function[] = noisv11(D,IT)

% Noisv11 calculate the total noise for a wire dipole antenna.
% [] = noisv11 (D, IT) Calculations is based on the distance to the
% sun and the wanted biased current and use the function ion, electron, shot and data2. 
%
% The result of the calculation will be plotted in a logarithmic scale.
%
% see also shot, ion, electron.
%
% $Id: nose.m,v 1.1 2013/04/10 14:52:00 Pansar Exp $

%% Declaration of constants
Units=irf_units;
Me= Units.me;      %% Electron mass            (Kg)
Mp= Units.mp;      %% Proton mass              (Kg)
eps0= Units.eps0;  %% Electric constant        (F/m)
qe= Units.e;       %% Elemetary charge         (C)
KB= Units.kB;      %% Boltsman konstatn        (J/K)
L=5;               %% Antenna lengths          (m)
r=0.575e-2;        %% Antenna radiens          (m)

tic;
Vsweep=[-5:0.00001:15];
P=[50 95 95];
%% Data gathering 
[Ne,Tp,B,V,RBt,RBe,RBi,RBp,RUt,RUe,RUi,RUp,VU,VB,Is,Ib,Iu,IT0]...
    = wp.data2(D,IT,P,Vsweep,r);

%%Tp  Proton Temperatur                     	(K)
%%RBe Antenna resistance biased electron     	(ohm)
%%RBi Antenna resistance biased proton       	(ohm)        
%%RBp Antenna resistance biased potonelectron	(ohm)
%%RBt Antenna resistance biased total        	(ohm) 
%%RUe Antenna resistance unbiased electron   	(ohm)
%%RUi Antenna resistance unbiased proton     	(ohm)        
%%RUp Antenna resistance unbiased potonelectron (ohm)
%%RUt Antenna resistance unbiased total      	(ohm)
%%UV Unbiased voltage                           (V)
%%BV Biased voltage                         	(V)
%%Is current saturation                     	(A)
%%Ib current biased                             (A)
%%Iu current biased                             (A)

Ne=Ne.*1e6;         %% Electron density             (m^-3) 
Te=Tp./2.5;         %% Electron temperatur          (K)
V=V.*1e3;           %% Solar wind velocity          (m/s)
B=B.*1e-9;          %% Magnetic field               (T)

Cp=30e-12;          %% Plasma capacitance          (F)

%% Calculations
for i=1:length(P),
%% Calculations for Issad
Fp(i)=sqrt(Ne(i)*qe.^2/(Me*eps0))/(2*pi);   %% Plasma frequency             (Hz)
LFe(i)=qe*B(i)/(2*pi*Me);                   %% Larmor frequency electron    (Hz)
LFp(i)=qe*B(i)/(2*pi*Mp);                   %% Larmor frequency proton      (Hz)
v(i)=sqrt((KB*Te(i))/Me);                   %% Electron thermal velocity    (m/s)
Ld(i)=sqrt(eps0*KB*Te(i)/(Ne(i)*(qe)^2));   %% Debye length                 (m)


f{i}=10.^[log10(10^-1):0.01:log10(1e6)];            % Frequency range

I{i}=wp.C.ion(f{i},Ne(i),Te(i),Tp(i),V(i),L);       % Ion QTN

E{i}=wp.C.electron(f{i},Ne(i),Te(i),L);             % Electron QTN  

SBe{i}=wp.C.shote(f{i},Ne(i),Te(i),RBt(i),Cp,L,Me); % Shot noise biased

SUe{i}=wp.C.shote(f{i},Ne(i),Te(i),RUt(i),Cp,L,Me); % Shot noise unbiased

n=1;
for j=f{i}
    if j<=LFe(i)
        E{i}(n)=0;
        n=n+1;
    end
end % Removing electron QTN values lower than the electron gyrofrequency

n=1;
for j=f{i}
    if j<=LFp(i)
        I{i}(n)=0;
        n=n+1;
    end
end  % Removing ion QTN value lower than the ion gyrofrequency

     
T{i}=sqrt(I{i}.^2+E{i}.^2); % Total QTN


%% Calculating the other noises 
IB(i)=sqrt(Is(i)^2+(Is(i)-Ib(i))^2); 
IU(i)=sqrt(Is(i)^2+(Is(i)-Iu(i))^2);

[Veb{i},Vib{i},Vpb{i},Vb{i},VEb{i},Vab{i}]...
    = wp.R.noisR(f{i},Te(i),Tp(i),RBe(i),RBi(i),RBp(i),IB(i),L,T{i});

[Veu{i},Viu{i},Vpu{i},Vu{i},VEu{i},Vau{i}] ...
    = wp.R.noisR(f{i},Te(i),Tp(i),RUe(i),RUi(i),RUp(i),IU(i),L,T{i});

 Un{i}=sqrt(VEu{i}.^2+Vu{i}.^2+Vau{i}.^2+SUe{i}.^2);
 Bi{i}=sqrt(VEb{i}.^2+Vb{i}.^2+Vab{i}.^2+SBe{i}.^2);
end

for i=1:length(P),
%% Plots
% Ploting the QTN
issad=figure(length(P)+i);
%Setup of plot
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

loglog(f{i},I{i},'-r',f{i},E{i},'-g',f{i},T{i},'-.k')

xlabel('Frequency (Hz)')
ylabel('VPS ((V/m)^2/Hz)')

line ([LFp(i) LFp(i)],[1e-20 1e-12],'color','red','LineStyle',':','LineWidth',2)
line ([LFe(i) LFe(i)],[1e-20 1e-12],'color','green','LineStyle',':','LineWidth',2)
line ([Fp(i) Fp(i)],[1e-20 1e-12],'color','cyan','LineStyle',':','LineWidth',2)

%%set(gca,'ylim',[10^-19 10^-13]);

legend('Ion ','Electron','Total','Larmor frequency ion',...
    'Larmor frequency electron','Plasma frequency',...
    'Location','Best')
grid on
xlim([10^-1 10^6])
ylim([10^-19 10^-11])
set(issad,'color','white'); % White background for the figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\issad',num2str(D*100,'%6.4g'),num2str(i... 
    ,'%6.4g'),'.eps']; %Location of saved plot
print( '-depsc2' , name ) % Saving plot


%% Plot Unbiased
unbias=figure(2*length(P)+i);
%Setup of plot
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
set(gcf,'paperpositionmode','auto') % To get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

loglog( f{i},Viu{i},'-r',f{i},Veu{i},'-g',f{i},Vpu{i},'-b',f{i},SUe{i},'-y',...
    f{i},Vau{i},'-m',f{i},VEu{i},'-c',...
    f{i},Un{i},'-.k')

xlabel('Frequency (Hz)')
ylabel('(V/m)^2/Hz')

legend('Thermal noise (Ion) ','Thermal noise (Electron) ','Thermal noise (Photoelectron)','Shot noise ',...
    'Amplifier noise','QTN','Total noise',...
    'Location','Best')

G=round(0.5+log10(sqrt(VEu{i}(1).^2+Vu{i}(1).^2+Vau{i}(1).^2)));
grid on
xlim([10^2 10^6])
ylim([10^-19 10^-11])
set(unbias,'color','white'); % White background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Unbiased',num2str(D*100,'%6.4g'),... 
    num2str(i,'%6.4g'),'.eps']; %Location of saved plot
print( '-depsc2' , name ) % Saving plot

%% Plot biased
bias=figure(3*length(P)+i);
%Setup of plot
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

loglog( f{i},Vib{i},'-r',f{i},Veb{i},'-g',f{i},Vpb{i},'-b',f{i},SBe{i},'-y',....
    f{i},Vab{i},'-m',f{i},VEb{i},'-c',f{i},Bi{i},'-.k')


xlabel('Frequency (Hz)')
ylabel('(V/m)^2/Hz')


legend('Thermal noise (Ion) ','Thermal noise (Electron) ','Thermal noise (Photoelectron)','Shot noise ',...
    'Amplifiern noise','QTN','Total noise',...
    'Location','Best')
G=round(0.5+log10(sqrt(VEb{i}(1).^2+Vb{i}(1).^2+Vab{i}(1).^2)));
grid on
xlim([10^2 10^6])
ylim([10^-19 10^-11])
set(bias,'color','white'); % White background for figures (default is grey)
name=['C:\Users\Modesty\Desktop\pic\Biased',num2str(D*100,'%6.4g'),... 
    num2str(i,'%6.4g'),'.eps']; %Location of saved plot
print( '-depsc2' , name ) % Saving plot

clf(figure(4*length(P)+i));
%% valuse from the calulations
num=figure(4*length(P)+i);
%Setup of plot
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
set(gcf,'paperpositionmode','auto') % To get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
set(gcf,'paperpositionmode','auto')

txstr(1)={['Distance to the sun =',num2str(D,'%6.4g'),'(AU)']};

txstr(2)={['Ne = ',num2str(Ne(i),'%6.4g'),' (m^3)']};
txstr(3)={['Te = ' ,num2str(8.61734e-5*Te(i),'%6.4g'),' (K)']};
txstr(4)={['Tp = ',num2str(8.61734e-5*Tp(i),'%6.4g'),' (K)']};
txstr(5)={['V = ',num2str(V(i),'%6.4g'),' (m/s)']};
txstr(6)={['B = ',num2str(B(i),'%6.4g'),' (T)']};
txstr(7)={['Kappa = 2']};

txstr(8)={['R total plasma biased = ',num2str(RBt(i),'%6.4g'),' (ohm)']};
txstr(9)={['R total plasma unbiased = ',num2str(RUt(i),'%6.4g'),' (ohm)']};
txstr(10)={['R ion plasma biased = ',num2str(RBi(i),'%6.4g'),' (ohm)']};
txstr(11)={['R ion plasma unbiased = ',num2str(RUi(i),'%6.4g'),' (ohm)']};
txstr(12)={['R electrons plasma biased = ',num2str(RBe(i),'%6.4g'),' (ohm)']};
txstr(13)={['R electrons plasma unbiased = ',num2str(RUe(i),'%6.4g'),' (ohm)']};
txstr(14)={['R photoelectron biased = ',num2str(RBp(i),'%6.4g'),' (ohm)']};
txstr(15)={['R photoelectron unbiased = ',num2str(RUp(i),'%6.4g'),' (ohm)']};
txstr(16)={['Voltage when unbiased = ',num2str(VU(i),'%6.4g'),' (V)']};
txstr(17)={['Voltage when biased = ',num2str(VB(i),'%6.4g'),' (V)']};
txstr(18)={['Current when unbiased = ',num2str(Iu(i),'%6.4g'),' (A)']};
txstr(19)={['Target bias current ',num2str(IT,'%6.4g'),' closest value ',...
    num2str(Ib(i),'%6.4g'),' (A)']};

txstr(20)={['Larmor frequency electron = ',num2str(LFe(i),'%6.4g'),' (Hz)']};
txstr(21)={['Larmor frequency ion = ',num2str(LFp(i),'%6.4g'),' (Hz)']};
txstr(22)={['Plasma frequency = ',num2str(Fp(i),'%6.4g'),' (Hz)']};
txstr(23)={['Debay length = ',num2str(Ld(i),'%6.4g'),' (m)']};

txstr(24)={['Frequency = ',num2str(f{i}(1),'%6.4g'),' to',num2str(f{i}(length(f{i}))...
    ,'%6.4g'),' (Hz)']};
txstr(25)={['Antenna length = ',num2str(L,'%6.4g'),'(m)']};

text(0.1e0,0.5e0,txstr,'HorizontalAlignment','left')
set(num,'color','white'); % White background for figures (default is grey)
set(gcf,'paperpositionmode','auto') % To get the same on paper as on screen
% to get bitmap file
name=['C:\Users\Modesty\Desktop\pic\data',num2str(D*100,'%6.4g'),... 
   num2str(i,'%6.4g'),'.png'];% Location of saved plot
print( '-dpng' , name ) % Saving plot
end
 toc   
end


