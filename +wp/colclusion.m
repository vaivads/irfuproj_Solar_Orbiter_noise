function[] = colclusion()
%%low
Fi=[0.01 100 10000];
I=[10^-12 10^-16 10^-16];

F=[0.100000000000000,1,10,100,1000,10000,100000,1000000];
B3=[2.83717387168589e-15,3.01502216931769e-15,3.01502164262486e-15,3.01496897431811e-15,3.00971189066986e-15,2.57538221726621e-15,4.81256618575455e-16,2.77760247412303e-18];
U3=[2.49625581476963e-13,2.49627176588159e-13,2.49578234119818e-13,2.44821869862773e-13,1.06355442106766e-13,2.48641250277699e-15,4.57899297900630e-16,2.27244668492770e-18];

figure(1)
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

loglog(Fi,I,'R-',F,U3,'B-.',F,B3,'B-')

grid off
xlim([10^-1 10^4])
ylim([10^-19 10^-3])
str={'Shocks; reconnection;','ion acoustic, ion cyclontron','and magnetosonic waves;','whistler turbulence'};
annotation('textbox',[0.28, 0.3, 0.55,0.53 ],'string',str,'FontSize',25,'HorizontalAlignment', 'center')

%%high
Fi=[1e3 1e4 5.5e5 3e6];
I=[10^-16 10^-18 10^-19 5e-19];

figure(2)
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

loglog(Fi,I,'R-',F,U3,'B-.',F,B3,'B-')

grid off
xlim([10^3 10^6])
ylim([10^-20 10^-8])
str1={'Thermal noise'};
annotation('textbox',[0.135, 0.25, 0.73,0.4 ],'string',str1,'FontSize',25)

str2={'Langmuir waves'};
annotation('textbox',[0.42, 0.4, 0.3,0.4 ],'string',str2,'FontSize',25,'HorizontalAlignment', 'center')
