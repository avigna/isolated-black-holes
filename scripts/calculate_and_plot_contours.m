function calculate_and_plot_contours(filename, metallicity, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
label = num2str(metallicity);
% MACROS
% Stellar types as used in COMPAS, originally defined by Hurley+2000
MASSIVE_MAIN_SEQUENCE   = 1;
NEUTRON_STAR            = 13;
BLACK_HOLE              = 14;
AgeOfUniverseMyr        = 13.772*1000;
RsolToAU                = 0.00465047;

% Data observed
% Sahu: https://ui.adsabs.harvard.edu/abs/2022arXiv220113296S/abstract
% Mass: 7.1 +/- 1.3 M_Sun
% transverse space velocity of ~45 km/s
mass_Sahu = 7.1;
error_Sahu = 1.3;
velocity_Sahu = 45;

% Lam: https://arxiv.org/abs/2202.01903
% Mass: 1.6 - 4.2 Msun
% Transverse motion: <25 km/s
mass_Lam = 2.9;
error_Lam = 1.3;
velocity_Lam = 25;

% Get isolated systems by resampling the primary of wide COMPAS binaries
[  sys_realMassZAMS, ...
   iso_real_total_mass_ZAMS,...
   iso_kick_magnitude,...
   iso_mass_SN,...
   iso_stellar_type_SN] = resampling_COMPAS(filename,debugFlag,lowerLimit);
realMassZAMS = sys_realMassZAMS+iso_real_total_mass_ZAMS;

% /BSE_Supernovae
SNe_component_speed_SN      = h5read(filename,'/BSE_Supernovae/ComponentSpeed(SN)');
SNe_eccentricity            = h5read(filename,'/BSE_Supernovae/Eccentricity');
SNe_mass_CP                 = h5read(filename,'/BSE_Supernovae/Mass(CP)');
SNe_mass_SN                 = h5read(filename,'/BSE_Supernovae/Mass(SN)');
SNe_SemiMajorAxis           = h5read(filename,'/BSE_Supernovae/SemiMajorAxis').*RsolToAU;
SNe_stellar_type_CP         = h5read(filename,'/BSE_Supernovae/Stellar_Type(CP)');
SNe_stellar_type_SN         = h5read(filename,'/BSE_Supernovae/Stellar_Type(SN)');
SNe_supernova_state         = h5read(filename,'/BSE_Supernovae/Supernova_State');
SNe_systemic_speed          = h5read(filename,'/BSE_Supernovae/SystemicSpeed');
SNe_unbound                 = h5read(filename,'/BSE_Supernovae/Unbound');

% SNe_supernova_state
% No supernova = 0
% Star 1 is the supernova = 1
% Star 2 is the supernova = 2
% Both stars are supernovae = 3

% Calculate quantities
SNe_total_mass      = SNe_mass_CP + SNe_mass_SN;
SNe_periapsis       = SNe_SemiMajorAxis.*(1-SNe_eccentricity);

% Values for creatings contours
% [bandwith, density, X, Y] = kde2d(data);
% Template
% data = [SNe_mass_SN(index_), SNe_component_speed_SN(index_)];
% [bandwith, density, X, Y] = kde2d(data,n,MIN_XY,MAX_XY);
% v =max(max(density)).*[0.01,1.0];
n=2^10;
MIN_XY = [1, 0];
MAX_XY = [100, 1000];

% Find indices of interest
% 0. Isolated BH
% 1. First born BH (when disrupted),
% 2. first born NS (when disrupted),
% 3. second born BH (when disrupted),
% 4. second born NS (when disrupted),
% 5. first born BH, when not disrupted (X-ray binarieS at birth),
% 6. first born NS, when not disrupted (X-ray binarieS at birth),
% 7. formation of bound NS-BH or BH-NS,
% 8. formation of bound BBH,
% 9. formation of bound NSNS,
% 10. NS-BH, BH-NS merger product,
% 11. BBH merger product,
% 12. NSNS merger product.

% 0. Isolated BH
index_zero      = find(iso_stellar_type_SN == BLACK_HOLE);
number_zero     = length(index_zero);
data_zero       = [iso_mass_SN(index_zero), iso_kick_magnitude(index_zero)];
[bandwith_zero, density_zero, X_zero, Y_zero] = kde2d(data_zero);
v_zero          = max(max(density_zero)).*[1.0-densityInverval,1.0];


% 1. First born BH (when disrupted)
index_one       = find( SNe_stellar_type_SN == BLACK_HOLE & ...
                        (SNe_supernova_state == 1 | SNe_supernova_state == 3) & ...
                        SNe_unbound == 1);
number_one      = length(index_one);
data_one        = [SNe_mass_SN(index_one), SNe_component_speed_SN(index_one)];
[bandwith_one, density_one, X_one, Y_one] = kde2d(data_one,n,MIN_XY,MAX_XY);
v_one           = max(max(density_one)).*[1.0-densityInverval,1.0];


% 2. First born NS (when disrupted)
index_two       = find( SNe_stellar_type_SN == NEUTRON_STAR & ...
                        (SNe_supernova_state == 1 | SNe_supernova_state == 3) & ...
                        SNe_unbound == 1);
number_two      = length(index_two);
data_two        = [SNe_mass_SN(index_two), SNe_component_speed_SN(index_two)];
[bandwith_two, density_two, X_two, Y_two] = kde2d(data_two,2^4,MIN_XY,[3 1000]);
v_two           = max(max(density_two)).*[1.0-densityInverval,1.0];

% 3. second born BH (when disrupted)
index_three     = find( (SNe_stellar_type_SN == BLACK_HOLE & SNe_supernova_state == 2 & SNe_unbound == 1) | ...
                        (SNe_stellar_type_CP == BLACK_HOLE & SNe_supernova_state == 3 & SNe_unbound == 1));
number_three    = length(index_three);
data_three      = [SNe_mass_SN(index_three), SNe_component_speed_SN(index_three)];
[bandwith_three, density_three, X_three, Y_three] = kde2d(data_three,n,MIN_XY,MAX_XY);
v_three         = max(max(density_three)).*[1.0-densityInverval,1.0];
 
% 5. first born BH, when not disrupted (X-ray binarieS at birth)
index_five      = find( SNe_stellar_type_SN == BLACK_HOLE & ...
                        SNe_stellar_type_CP == MASSIVE_MAIN_SEQUENCE & ...
                        SNe_unbound == 0);
number_five     = length(index_five);
data_five       = [SNe_total_mass(index_five), SNe_systemic_speed(index_five)];
[bandwith_five, density_five, X_five, Y_five] = kde2d(data_five);
v_five          = max(max(density_five)).*[1.0-densityInverval,1.0];

% 6. first born NS, when not disrupted (X-ray binarieS at birth)
index_six       = find( SNe_stellar_type_SN == NEUTRON_STAR &...
                        SNe_stellar_type_CP == MASSIVE_MAIN_SEQUENCE &...
                        SNe_unbound == 0);
number_six      = length(index_six);
data_six = [SNe_total_mass(index_six), SNe_systemic_speed(index_six)];
[bandwith_six, density_six, X_six, Y_six] = kde2d(data_six,2^9,[1.5 0],[107 500]);
v_six           = max(max(density_six)).*[1.0-densityInverval,1.0];

% 10. NS-BH, BH-NS mergers
index_all_BHNSs = find( ((SNe_stellar_type_SN == NEUTRON_STAR & SNe_stellar_type_CP == BLACK_HOLE) |...
                        (SNe_stellar_type_SN == BLACK_HOLE & SNe_stellar_type_CP == NEUTRON_STAR))...
                        & SNe_unbound == 0);
if debugFlag
    length(index_all_BHNSs)
end
SNe_total_mass_ten_temp         = SNe_total_mass(index_all_BHNSs);
SNe_systemic_speed_ten_temp     = SNe_systemic_speed(index_all_BHNSs);
SNe_merger_time_GWs_ten_temp    = calculateMergerTimescale(SNe_mass_SN(index_all_BHNSs),SNe_mass_CP(index_all_BHNSs),SNe_SemiMajorAxis(index_all_BHNSs),SNe_eccentricity(index_all_BHNSs),0);
index_ten       = find(SNe_merger_time_GWs_ten_temp<AgeOfUniverseMyr);
number_ten      = length(index_ten);
SNe_total_mass_ten      = SNe_total_mass_ten_temp(index_ten);
SNe_systemic_speed_ten  = SNe_systemic_speed_ten_temp(index_ten);

data_ten        = [SNe_total_mass_ten, SNe_systemic_speed_ten];
[bandwith_ten, density_ten, X_ten, Y_ten] = kde2d(data_ten,n,MIN_XY,MAX_XY);
v_ten           = max(max(density_ten)).*[1.0-densityInverval,1.0];

SNe_SemiMajorAxis_ten_temp  = SNe_SemiMajorAxis(index_all_BHNSs);
SNe_SemiMajorAxis_ten       = SNe_SemiMajorAxis_ten_temp(index_ten);
% SNe_periapsis_ten_temp      = SNe_periapsis(index_all_BHNSs);
% SNe_periapsis_ten           = SNe_periapsis_ten_temp(index_ten);

% 11. BH-BH mergers
index_all_BBHs     = find( ((SNe_stellar_type_SN == BLACK_HOLE & SNe_stellar_type_CP == BLACK_HOLE))...
                        & SNe_unbound == 0);
if debugFlag
    length(index_all_BBHs)
end
SNe_total_mass_eleven_temp         = SNe_total_mass(index_all_BBHs);
SNe_systemic_speed_eleven_temp     = SNe_systemic_speed(index_all_BBHs);
SNe_merger_time_GWs_eleven_temp    = calculateMergerTimescale(SNe_mass_SN(index_all_BBHs),SNe_mass_CP(index_all_BBHs),SNe_SemiMajorAxis(index_all_BBHs),SNe_eccentricity(index_all_BBHs),0);
index_eleven       = find(SNe_merger_time_GWs_eleven_temp<AgeOfUniverseMyr);
number_eleven      = length(index_eleven);
SNe_total_mass_eleven = SNe_total_mass_eleven_temp(index_eleven);
SNe_systemic_speed_eleven = SNe_systemic_speed_eleven_temp(index_eleven);

data_eleven        = [SNe_total_mass_eleven, SNe_systemic_speed_eleven];
[bandwith_eleven, density_eleven, X_eleven, Y_eleven] = kde2d(data_eleven,2^8,MIN_XY,MAX_XY);
v_eleven           = max(max(density_eleven)).*[1.0-densityInverval,1.0];

SNe_SemiMajorAxis_eleven_temp       = SNe_SemiMajorAxis(index_all_BBHs);
SNe_SemiMajorAxis_eleven            = SNe_SemiMajorAxis_eleven_temp(index_eleven);
% SNe_periapsis_eleven_temp           = SNe_periapsis(index_all_BBHs);
% SNe_periapsis_eleven                = SNe_periapsis_eleven_temp(index_eleven);

% 12. NS-NS mergers
index_all_BNSs      = find( ((SNe_stellar_type_SN == NEUTRON_STAR & SNe_stellar_type_CP == NEUTRON_STAR))...
                        & SNe_unbound == 0);
if debugFlag
    length(index_all_BNSs)
end
SNe_total_mass_twelve_temp         = SNe_total_mass(index_all_BNSs);
SNe_systemic_speed_twelve_temp     = SNe_systemic_speed(index_all_BNSs);
SNe_merger_time_GWs_twelve_temp    = calculateMergerTimescale(SNe_mass_SN(index_all_BNSs),SNe_mass_CP(index_all_BNSs),SNe_SemiMajorAxis(index_all_BNSs),SNe_eccentricity(index_all_BNSs),0);
index_twelve        = find(SNe_merger_time_GWs_twelve_temp<AgeOfUniverseMyr);
number_twelve       = length(index_twelve);
SNe_total_mass_twelve = SNe_total_mass_twelve_temp(index_twelve);
SNe_systemic_speed_twelve = SNe_systemic_speed_twelve_temp(index_twelve);

data_twelve         = [SNe_total_mass_twelve, SNe_systemic_speed_twelve];
[bandwith_twelve, density_twelve, X_twelve, Y_twelve] = kde2d(data_twelve,2^6,MIN_XY,[5, 500]);
v_twelve            = max(max(density_twelve)).*[1.0-densityInverval,1.0];

SNe_SemiMajorAxis_twelve_temp       = SNe_SemiMajorAxis(index_all_BNSs);
SNe_SemiMajorAxis_twelve            = SNe_SemiMajorAxis_twelve_temp(index_twelve);
% SNe_periapsis_twelve_temp           = SNe_periapsis(index_all_BNSs);
% SNe_periapsis_twelve                = SNe_periapsis_twelve_temp(index_twelve);

if debugFlag
    number_zero
    number_one
    number_two
    number_three
    number_five
    number_six
    number_ten
    number_eleven
    number_twelve
end

% Make histogram of separation
hist                    = histogram(log10(SNe_SemiMajorAxis(index_five)),'Normalization','pdf');
bins_separation_five    = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
values_separation_five  = hist.Values;

hist                    = histogram(log10(SNe_SemiMajorAxis(index_six)),'Normalization','pdf');
bins_separation_six     = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
values_separation_six   = hist.Values;

hist                    = histogram(log10(SNe_SemiMajorAxis_ten),'Normalization','pdf');
bins_separation_ten     = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
values_separation_ten   = hist.Values;

hist                    = histogram(log10(SNe_SemiMajorAxis_eleven),'Normalization','pdf');
bins_separation_eleven  = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
values_separation_eleven= hist.Values;

hist                    = histogram(log10(SNe_SemiMajorAxis_twelve),'Normalization','pdf');
bins_separation_twelve  = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
values_separation_twelve= hist.Values;

% Will delete this if referee doesn't require it
% % Make histogram of periapsis
% hist                    = histogram(log10(SNe_periapsis(index_five)),'Normalization','pdf');
% bins_periapsis_five     = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
% values_periapsis_five   = hist.Values;
% 
% hist    = histogram(log10(SNe_periapsis(index_six)),'Normalization','pdf');
% bins_periapsis_six      = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
% values_periapsis_six    = hist.Values;
% 
% hist    = histogram(log10(SNe_periapsis_ten),'Normalization','pdf');
% bins_periapsis_ten      = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
% values_periapsis_ten    = hist.Values;
% 
% hist = histogram(log10(SNe_periapsis_eleven),'Normalization','pdf');
% bins_periapsis_eleven   = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
% values_periapsis_eleven = hist.Values;
% 
% hist = histogram(log10(SNe_periapsis_twelve),'Normalization','pdf');
% bins_periapsis_twelve   = hist.BinEdges(1:end-1)+(hist.BinWidth/2);
% values_periapsis_twelve = hist.Values;

clear hist
clf
% Plot
if plotFlag
    figure()
    fs=18;
    lw=2.0;
    sz=1.0;
    solar=char(9737);
    string1=['Mass [M_',solar,']'];  
    string2='CoM Speed [km s^{-1}]';
    alphaNum = 0.4;
    Ylim = [0 400];
    Xlim = [1 100];

    color0 = 105.*[1 1 1]./255;
    color1 = [0    0.4470    0.7410];
    color2 = [    0.8500    0.3250    0.0980];
    color3 = [    0.9290    0.6940    0.1250];
    color4 = [    0.4940    0.1840    0.5560];
    color5 = [    0.4660    0.6740    0.1880];
    color6 = [    0.3010    0.7450    0.9330];
    color7 = [    0.6350    0.0780    0.1840];

    % Figure 2a
    hold on
    text(1.15,375,strcat('Z=',label),'FontSize',fs,'FontName','Times New Roman')
    ylim(Ylim)
    xlim(Xlim)
    set(gca,'XScale','log')
    xticks([1 10 100])
    xticklabels({'1','10','100'})
    xlabel(string1,'FontSize',fs)
    ylabel(string2,'FontSize',fs)
    ax=gca;
    ax.FontSize=fs;
    ax.FontName = 'Times New Roman';
    box on

    contour(X_zero,Y_zero,density_zero,v_zero,'Color',color0,'Linewidth',lw)
    contour(X_one,Y_one,density_one,v_one,'Color',color1,'Linewidth',lw)
    contour(X_three,Y_three,density_three,v_three,'Color',color2,'Linewidth',lw)
    contour(X_five,Y_five,density_five,v_five,'Color',color3,'Linewidth',lw)
    contour(X_six,Y_six,density_six,v_six,'Color',color4,'Linewidth',lw)

    errorbar(mass_Sahu,velocity_Sahu,error_Sahu,'horizontal','k','LineWidth',lw,'HandleVisibility','off')
    errorbar(mass_Lam,velocity_Lam,error_Lam,'horizontal','k','LineWidth',lw,'HandleVisibility','off')
    annotation('textarrow',[0.73 0.46],[0.5 0.24],'String','Sahu et al.','Fontsize',fs,'FontName','Times New Roman')
    annotation('textarrow',[0.21 0.27],[0.5 0.20],'String','Lam et al.','Fontsize',fs,'FontName','Times New Roman')
    

    if debugFlag
        scatter(iso_mass_SN(index_zero), iso_kick_magnitude(index_zero), sz, color0, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')        
        scatter(SNe_mass_SN(index_one), SNe_component_speed_SN(index_one), sz, color1, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
        scatter(SNe_mass_SN(index_three), SNe_component_speed_SN(index_three), sz, color2, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
        scatter(SNe_total_mass(index_five), SNe_systemic_speed(index_five), sz, color3, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
        scatter(SNe_total_mass(index_six), SNe_systemic_speed(index_six), sz, color4, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
    end

    legend( 'BH_{sin}',...
            'BH_1',...
            'BH_2',...
            'BH-MS',...
            'NS-MS',...
            'box','off',...
            'location','northeast',...
            'FontName','Times New Roman')
    if savePlotFlag
        print(gcf,strcat('../plots/png/figure2a_Z=',label,'.png'),'-dpng','-r300');
        saveas(gcf,strcat('../plots/pdf/figure2a_Z=',label,'.pdf'))
    end

    % Figure 2b
    figure()
    hold on
    text(1.15,375,strcat('Z=',label),'FontSize',fs,'FontName','Times New Roman')
    ylim(Ylim)
    xlim(Xlim)
    set(gca,'XScale','log')
    xticks([1 10 100])
    xticklabels({'1','10','100'})
    xlabel(string1,'FontSize',fs)
    ylabel(string2,'FontSize',fs)
    ax=gca;
    ax.FontSize=fs;
    ax.FontName = 'Times New Roman';    
    box on

    contour(X_ten,Y_ten,density_ten,v_ten,'Color',color5,'Linewidth',lw)
    contour(X_eleven,Y_eleven,density_eleven,v_eleven,'Color',color6,'Linewidth',lw)    
    contour(X_twelve,Y_twelve,density_twelve,v_twelve,'Color', color7,'Linewidth',lw)    

    errorbar(mass_Sahu,velocity_Sahu,error_Sahu,'horizontal','k','LineWidth',lw,'HandleVisibility','off')
    errorbar(mass_Lam,velocity_Lam,error_Lam,'horizontal','k','LineWidth',lw,'HandleVisibility','off')
    annotation('textarrow',[0.73 0.46],[0.5 0.24],'String','Sahu et al.','Fontsize',fs,'FontName','Times New Roman')
    annotation('textarrow',[0.21 0.27],[0.5 0.20],'String','Lam et al.','Fontsize',fs,'FontName','Times New Roman')


    if debugFlag    
        scatter(SNe_total_mass_ten, SNe_systemic_speed_ten, sz, color5, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
        scatter(SNe_total_mass_eleven, SNe_systemic_speed_eleven, sz, color6, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
        scatter(SNe_total_mass_twelve, SNe_systemic_speed_twelve, sz, color7, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')
    end

    legend( 'BH-NS',...
        'BH-BH',...
        'NS-NS',...
        'box','off',...
        'location','northeast')

    if savePlotFlag
        print(gcf,strcat('../plots/png/figure2b_Z=',label,'.png'),'-dpng','-r300');
        saveas(gcf,strcat('../plots/pdf/figure2b_Z=',label,'.pdf'))
    end

    if metallicity == 0.02
        % Figure 4    
        figure()
        lw2=2.5;
        hold on
        text(2.75,2.75,strcat('Z=',label),'FontSize',fs,'FontName','Times New Roman')
        xlim([-4 4])
        xticks([-4:4])
        xlabel('log_{10} (a [AU])','FontSize',fs,'FontName','Times New Roman')
        ylabel('Normalized distribution','FontSize',fs,'FontName','Times New Roman')
        ax=gca;
        ax.FontSize=fs;
        ax.FontName = 'Times New Roman';
        box on
    
        N=10;
        maxVal=3;
        alphaNum=0.1;
        x = linspace(log10(0.18),log10(230),N);
        y = repmat([0;maxVal],1,N);
        h1=fill([x flip(x)], [y(1,:) flip(y(2,:))],'k','EdgeColor','none','HandleVisibility','off');
        set(h1,'FaceAlpha',alphaNum)% if size(y0,1)==2 %plot shaded area       
        text(-0.5,2,{'Excluded region','Sahu et al.'},'Color','k','Fontsize',fs,'FontName','Times New Roman')  
    
        plot(bins_separation_five,values_separation_five,'-','Color',color3,'LineWidth',lw2)
        plot(bins_separation_six,values_separation_six,'-','Color',color4,'LineWidth',lw2)
        plot(bins_separation_ten,values_separation_ten,':','Color',color5,'LineWidth',lw2)
        plot(bins_separation_eleven,values_separation_eleven,':','Color',color6,'LineWidth',lw2)
        plot(bins_separation_twelve,values_separation_twelve,':','Color',color7,'LineWidth',lw2)
    
        legend( 'BH-MS', ...
                'NS-MS',...
                'BH-NS',...
                'BH-BH',...
                'NS-NS',...
                'box','off',...
                'location','northwest')
    
        if savePlotFlag
            print(gcf,strcat('../plots/png/figure7_Z=',label,'.png'),'-dpng','-r300');
            saveas(gcf,strcat('../plots/pdf/figure7_Z=',label,'.pdf'))
        end    
    end

end

% Gather and save data
M = [   realMassZAMS, ...
        number_zero, ...
        number_one, ...
        number_three, ...
        number_five, ...
        number_six, ...
        number_ten, ...
        number_eleven, ...
        number_twelve];

if savaDataFlag
    save(strcat('../data/Z_',label,'.mat'),'M')
end

if debugFlag
    test_plot_KDE(data_three, label, densityInverval)
end

end