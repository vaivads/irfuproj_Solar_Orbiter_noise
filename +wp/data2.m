function[Ne,Tp,Be,Ve,BRT,BRE,BRI,BRP,URT,URE,URI,URP,VU,VB,Is,Ib,Iu,IT0] = data2(D,IB,P,Vsweep,r)

% data2 calculate the plasma parameter as function of the distance to sun.
% The input data is collected from the Helios 1 and Helios 2 satellite.   
%
% [Ne,Te,Be,Ve,BRT,BRE,BRI,BRP,URT,URE,URI,URP,UV,BV,Is,Ib,Iu] = data2(D,IB,P,Vsweep,r)
% input:
% D is distances to sun .
% IB is the wanted current biased for the antenna.
% P is the average probability.
% Vsweep is the voltage sweep band
% r is the antenna radius.
%
% Output:
%Ne partial density. 
%Tp proton temperature. 
%Be magnetic field.
%Ve solar wind speed.
%BRT total resistance biased.
%BRE electron resistance biased.
%BRI Ion resistance biased.
%BRP photoelectron resistance biased.
%URT total resistance unbiased.
%URE electron resistance unbiased.
%URI Ion resistance unbiased.
%URP photoelectron resistance unbiased.
%VU  voltage when unbiased.
%VB  voltage when biased.
%Is  current when saturated.
%Ib  current when biased.
%Iu  current when unbiased.
%
% 
% see also: 
%
% $Id: nose.m,v 1.0 2013/04/10 14:48:00 Pansar Exp $
%
load 'wp/Helios_data'

i=1;
n=1;
x=helios1.heliocentricDistance(:,2);

while n<=length(x)
    
    if ((D-0.01)<x(n)) & ((D+0.01)>x(n))
          
        N(i)=helios1.protonDensity(n,2);
     
        T(i)=helios1.protonTemp(n,2);
       
        B(i)=helios1.B(n,2);
       
        V(i)=helios1.flowSpeed(n,2);
            
        i=i+1;
        
    end
    
    n=n+1;
end     % Collect values from Helios 1
n=1;
i=1;
x=helios2.heliocentricDistance(:,2);

while n<=length(x)   % Collect values from Helios 2
    
    if ((D-0.01)<x(n)) & ((D+0.01)>x(n))
          
        N(i)=helios2.protonDensity(n,2);
     
        T(i)=helios2.protonTemp(n,2);
       
        B(i)=helios2.B(n,2);
       
        V(i)=helios2.flowSpeed(n,2);
            
        i=i+1;
        
    end
    
    n=n+1;
end

for i = 1: length (P) % Removes NAN values
    
Ne(i) = prctile(N(~isnan(N)),P(i));
        
Tp(i) = prctile(T(~isnan(T)),P(i));
    
Be(i) = prctile(B(~isnan(B)),P(i));

Ve(i) = prctile(V(~isnan(V)),P(i));

if i==1  %Average case
  
elseif i==2 % High temperatur case
    n=1;
    c=1;
    L=[];
    K=[];
    M=[];  
    for j=T
       
       if Tp(i)-10000 <= j & j <= Tp(i)+10000
           L(c)=N(n);
           K(c)=B(n);
           M(c)=V(n);
           c=c+1;
       end
       n=n+1;
    end
   Ne(i) = mean(L(~isnan(L)));
   Be(i) = mean(K(~isnan(K)));
   Ve(i) = mean(M(~isnan(M)));
   
elseif i==3 % High density case
    n=1;
    c=1;
    L=[];
    K=[];
    M=[];
    for j=N
       
       if j <= Ne(i)+3 & j >= Ne(i)-3
           L(c)=T(n);
           K(c)=B(n);
           M(c)=V(n);
           c=c+1;
       end
       n=n+1;
    end
   Tp(i) = mean(L(~isnan(L)));
   Be(i) = mean(K(~isnan(K)));
   Ve(i) = mean(M(~isnan(M)));
   
end


   
end

for i=1:length(P),
    
    % Setup for probe_current
    probe.type  ='cylindrical';                             %'spherical','cylindrical','arbitrary'
    probe.surface ='gold';                                  %'themis','cassini' (one of flags in lp.photocurrent)
    probe.cross_section_area =r*2*5;                        %In m2
    probe.total_area =2*pi*(r*5);                           %In m2
    U_probe    =Vsweep;                                     %Probe potential (can be vector or matrix) 
    R_sun      =D;                                          %Distance from sun in AU
    UV_factor  =1;                                          %Default is 1
    %Describes plasma components (structure)
    plasma.q = [-1 1];                                      %- charge of species in e (the length of this vector corresponds to number of species)
    plasma.m = [0 1];                                       %Mass of species in proton masses (0 corresponds to e- mass)
    plasma.n = [Ne(i) Ne(i)];                               %Density of species [cc]
    plasma.T = [Tp(i)/2.5*8.61734e-5 Tp(i)*8.61734e-5];     %Temperature [eV]
    plasma.vsc =[Ve(i) Ve(i)];                              %Velocity of probe wrt. mmedia [m/s]
    
    [J_probe, J_photo, J_plasma]=lp.probe_current(probe,U_probe,R_sun,UV_factor,plasma);
    
    [V0 V0] = min(abs(Vsweep));
    
    Is(i)=min(J_photo);
    [idx idx] = min(abs(J_probe));

    if idx<21  % Check if the required values is on too high or too low
        idx=21;

        elseif idx>length(J_probe)-21
        idx=length(J_probe)-21;
    end

    
    % Calculating the resistance when it is unbiased 
    URT(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_probe(idx+20)-J_probe(idx-20));
    URP(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_photo(idx+20)-J_photo(idx-20));
    URE(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_plasma{1}(idx+20)-J_plasma{1}(idx-20));
    URI(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_plasma{2}(idx+20)-J_plasma{2}(idx-20));
    
    Iu(i)=J_probe(idx);
    VU(i)=U_probe(idx);
    IT0(i)=J_probe(V0);

    [idx idx] = min(abs(J_probe-IB));

    if idx<21
        idx=21;

    elseif idx>length(J_probe)-21
         idx=length(J_probe)-21;
    end
    % Calculating the resistance when it is biased 
    BRT(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_probe(idx+20)-J_probe(idx-20));
    BRP(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_photo(idx+20)-J_photo(idx-20));
    BRE(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_plasma{1}(idx+20)-J_plasma{1}(idx-20));
    BRI(i)=(U_probe(idx+20)-U_probe(idx-20))/(J_plasma{2}(idx+20)-J_plasma{2}(idx-20));   
    
    Ib(i)=J_probe(idx);
    VB(i)=U_probe(idx);
    
    % Plot of the collected values at  the distance
    PL=figure(1);
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

    plot(N,T,'+g',Ne,Tp,'+r')
    xlabel('Particle density  (m^3)')
    ylabel('Particle temperature (K)')
    legend('Data point','Calulation point','Location','Best')
    
    set(gcf,'color','white'); % White background for figures (default is grey)
    set(gcf,'paperpositionmode','auto') % To get the same on paper as on screen
    name=['C:\Users\Modesty\Desktop\pic\U-I\vs',num2str(D*100,'%6.4g'),...
        num2str(i,'%6.4g'),'.eps'];  % Location of saved plot
    print( '-depsc2' , name ) %  Saving plot
    
    % Plot of the U-I curve
    UI=figure((1+i));
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

    plot(U_probe,J_probe,'-k',U_probe,J_photo,'-b',U_probe,J_plasma{1},'-g',U_probe,J_plasma{2},'-r')
   
    ylim([-2.0e-5 0.5e-5])
    line ([VU(i) VU(i)],[0.5e-5 -2.0e-5],'color','magenta','LineStyle','-.','LineWidth',1.5)
    line ([VB(i) VB(i)],[0.5e-5 -2.0e-5],'color','cyan','LineStyle','-.','LineWidth',1.5)
    
    legend('Total','Photoelectrons','Electron','Ions','Unbiased','Biased','Location','Best')
    xlabel('Voltage (V)')
    ylabel('Current (A)')
    grid on
    set(gcf,'color','white'); % White background for figures (default is grey)
    set(gcf,'paperpositionmode','auto') % To get the same on paper as on screen
    name=['C:\Users\Modesty\Desktop\pic\U-I\UI',num2str(D*100,'%6.4g'),...
     num2str(i,'%6.4g'),'.eps']; % Location of saved plot
    print( '-depsc2' , name )% Saving plot
end

end





