function[] = conclud()
if 1, % calculate spectra 1AU
    Vf=630;T=60;B=15;nn=20; % R=1AU
    VA=20*B/sqrt(nn);
    SP.power(1)=2e5*1e-6;% power of signal at first frequency
    SP.Epower(1)=SP.power(1)*1e-6*VA^2;
    SP.EpowerVf(1)=SP.power(1)*1e-6*Vf^2;
    SP.R_AU=1; % distance in AU
    
    f_dop_ion=Vf/(2*pi*(100*sqrt(2*T)/B))*0.5;
    f_dop_electron=Vf/(2*pi*(sqrt(10*T)/B))*0.5;
    SP.R_RS=1/0.00465; % distance in AU
    SP.f=[f_range(1) 1e-3 f_dop_ion f_dop_electron f_dop_electron*10];
    SP.slopes=[-1 -1.6 -2.8 -4];
    for i=2:length(SP.f),
        SP.power(i)=SP.power(i-1)*10^(log10(SP.f(i)/SP.f(i-1))*SP.slopes(i-1));
        SP.Epower(i)=SP.power(i)*1e-6*VA^2;
        SP.EpowerVf(i)=SP.power(i)*1e-6*Vf^2;
    end
    SP1AU=SP;
end
if 1, % calculate spectra R=.3AU
    Vf=385;T=4;B=40;nn=25; % R=.3AU standard conditions
    VA=20*B/sqrt(nn);
    SP.power(1)=3e6*1e-6; % power of signal at first frequency
    SP.Epower(1)=SP.power(1)*1e-6*Vf^2;
    SP.R_AU=.3; % distance in AU
    f_dop_ion=Vf/(2*pi*(100*sqrt(2*T)/B))*0.5;
    f_dop_electron=Vf/(2*pi*(sqrt(10*T)/B))*0.5;
    SP.R_RS=1/0.00465; % distance in AU
    SP.f=[f_range(1) 6e-3 f_dop_ion f_dop_electron f_dop_electron*10];
    SP.slopes=[-1 -1.6 -2.8 -4];
    for i=2:length(SP.f),
        SP.power(i)=SP.power(i-1)*10^(log10(SP.f(i)/SP.f(i-1))*SP.slopes(i-1));
        SP.Epower(i)=SP.power(i)*1e-6*VA^2;
    end
    SP03AU=SP;
end
 
if 1, % initialize figure
    figure(11);clf
    h=irf_plot(1);
    set(h,'position',[0.15 0.1 0.75 0.75]);
    set(gcf,'defaultLineLineWidth',2);
    set(gcf,'PaperUnits','centimeters')
    xSize = 12; ySize = 13;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
end
if 1, % electric field plot
    hca=h(1);
    loglog(hca,SP1AU.f,SP1AU.Epower,'b.-','markersize',20);
    hold(hca,'on');
    loglog(hca,SP03AU.f,SP03AU.Epower,'r.-','markersize',20);
    
    set(hca,'xlim',f_range)
    set(hca,'ylim',PE_range)
    set(hca,'xtick',10.^[log10(f_range(1)):1:log10(f_range(2))]),
    set(hca,'ytick',10.^[log10(PE_range(1)):2:log10(PE_range(2))]),
    grid(hca,'on');
    set(hca,'xminorgrid','off');
    set(hca,'yminorgrid','off');
    
    ylabel(hca,'S_E [(V/m)^2/Hz]');
    xlabel(hca,'frequency [Hz]');
    
    text(0.97,0.8,'spectra at R=1 AU','fontsize',12,'fontweight','demi','color','b','units','normalized','horizontalalignment','right','parent',hca);
    text(0.97,0.85,'spectra at R=0.3 AU','fontsize',12,'fontweight','demi','color','r','units','normalized','horizontalalignment','right','parent',hca);
    
end