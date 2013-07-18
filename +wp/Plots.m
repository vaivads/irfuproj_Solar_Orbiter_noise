function[] = Plots()

load '/Users/wicpan/Dropbox/IRFU/matlab/helios/Helios_data'
x=helios1.heliocentricDistance(:,2);

     figure((1));
     plot(helios1.protonDensity(:,2),helios1.protonTemp(:,2),'+')
     xlabel('Partical densety(m^3)')
     ylabel('Temperatur (K)')

P=[50 95 95];     
p=1;
for D =[0.3 0.4 0.5 0.6 0.7 0.8 0.9]

N=0;
T=0;
B=0;
V=0;

i=1;
n=1;
while n<=length(x)
    
    if ((D-0.01)<x(n)) & ((D+0.01)>x(n))
          
        N(i)=helios1.protonDensity(n,2);
     
        T(i)=helios1.protonTemp(n,2);
       
        B(i)=helios1.B(n,2);
       
        V(i)=helios1.flowSpeed(n,2);
            
        i=i+1;
        
    end
    
    n=n+1;
end
n=1;
x=helios2.heliocentricDistance(:,2);

while n<=length(x)
    
    if ((D-0.01)<x(n)) & ((D+0.01)>x(n))
          
        N(i)=helios2.protonDensity(n,2);
     
        T(i)=helios2.protonTemp(n,2);
       
        B(i)=helios2.B(n,2);
       
        V(i)=helios2.flowSpeed(n,2);
            
        i=i+1;
        
    end
    
    n=n+1;
end


     

 for i = 1:length(P)
    
Ne(i) = prctile(N(~isnan(N)),P(i));
        
Tp(i) = prctile(T(~isnan(T)),P(i));
    
Be(i) = prctile(B(~isnan(B)),P(i));

Ve(i) = prctile(V(~isnan(V)),P(i));


if i==2
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

elseif i==3
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
     figure(1+p);
     set(0,'defaultLineLineWidth', 1.5);
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    set(gcf,'defaultAxesFontUnits','pixels');
    set(gcf,'defaultTextFontUnits','pixels');
    set(gcf,'PaperUnits','centimeters')
    xSize = 12; ySize = 12;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
    set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
    clear xLeft xSize sLeft ySize yTop
    set(gcf,'paperpositionmode','auto')
    set(gcf,'color','white'); % white background for figures (default is grey)
    
     plot(N,T,'+r',Ne,Tp,'+g')
     xlabel('Particle density(m^3)')
     ylabel('Particle temperature (K)')
     title([num2str(D*100 ,'%6.4g'), 'Au from the sun '])
     name=['\Users\wicpan\Dropbox\IRFU\pic\U-I\DT',num2str(D*100,'%6.4g'),'.eps'];
     print( '-depsc2' , name )
     
     figure(8+p);
     plot(B,T,'+b',Be,Tp,'+g')
     xlabel('Magnetic field')
     ylabel('Particletemperature (K)')
     title([num2str(D*100 ,'%6.4g'), 'Au from the sun '])

     
     figure(16+p);
     plot(V,T,'+k',Ve,Tp,'+g')
     xlabel('solar wind ')
     ylabel('Temperatur (K)')
     title([num2str(D*100 ,'%6.4g'), 'Au from the sun '])
    p=p+1; 
    N=[];
    T=[];
end