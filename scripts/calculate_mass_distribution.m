function calculate_mass_distribution(filename, metallicity, debugFlag, plotFlag, savePlotFlag, lowerLimit)

if debugFlag~=1
    debugFlag = 0;
end

if debugFlag
    dir
    debugFlag
end

% filename = strcat(dir,metallicityLabel,'/COMPAS_Output.h5');

% MACROS
% Stellar types as used in COMPAS, originally defined by Hurley+2000
NEUTRON_STAR            = 13;
BLACK_HOLE              = 14;

% Get isolated systems
[  sys_realMassZAMS, ...
   iso_real_total_mass_ZAMS,...
   iso_kick_magnitude,...
   iso_mass_SN,...
   iso_stellar_type_SN] = resampling_COMPAS(filename,debugFlag,lowerLimit);

% realMassZAMS = sys_realMassZAMS+iso_real_total_mass_ZAMS;

% /BSE_Supernovae
% SNe_                 = h5read(filename,'/BSE_Supernovae/');
SNe_component_speed_SN      = h5read(filename,'/BSE_Supernovae/ComponentSpeed(SN)');
SNe_mass_SN                 = h5read(filename,'/BSE_Supernovae/Mass(SN)');
SNe_stellar_type_CP         = h5read(filename,'/BSE_Supernovae/Stellar_Type(CP)');
SNe_stellar_type_SN         = h5read(filename,'/BSE_Supernovae/Stellar_Type(SN)');
SNe_supernova_state         = h5read(filename,'/BSE_Supernovae/Supernova_State');
SNe_unbound                 = h5read(filename,'/BSE_Supernovae/Unbound');

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

% No supernova = 0
% Star 1 is the supernova = 1
% Star 2 is the supernova = 2
% Both stars are supernovae = 3

% 0. Isolated BH
index_zero      = find(iso_stellar_type_SN == BLACK_HOLE);
number_zero     = length(index_zero);
data_zero       = [iso_mass_SN(index_zero), iso_kick_magnitude(index_zero)];

h_zero              = histogram(data_zero(:,1));
h_zero_values       = h_zero.Values;
h_zero_bin_center   = h_zero.BinEdges(1:end-1)+(h_zero.BinWidth/2);
h_zero_bin_edges    = h_zero.BinEdges;

% 00. Isolated NS
index_zerozero      = find(iso_stellar_type_SN == NEUTRON_STAR);
number_zerozero     = length(index_zerozero);
data_zerozero       = [iso_mass_SN(index_zerozero), iso_kick_magnitude(index_zerozero)];

% 1. First born BH (when disrupted)
index_one       = find( SNe_stellar_type_SN == BLACK_HOLE & ...
                        (SNe_supernova_state == 1 | SNe_supernova_state == 3) & ...
                        SNe_unbound == 1);
number_one      = length(index_one);
data_one        = [SNe_mass_SN(index_one), SNe_component_speed_SN(index_one)];

h_one              = histogram(data_one(:,1),h_zero_bin_edges);
h_one_values       = h_one.Values;
h_one_bin_center   = h_one.BinEdges(1:end-1)+(h_one.BinWidth/2);

% 2. First born NS (when disrupted)
index_two       = find( SNe_stellar_type_SN == NEUTRON_STAR & ...
                        (SNe_supernova_state == 1 | SNe_supernova_state == 3) & ...
                        SNe_unbound == 1);
number_two      = length(index_two);
data_two        = [SNe_mass_SN(index_two), SNe_component_speed_SN(index_two)];

% 3. second born BH (when disrupted)
index_three     = find( (SNe_stellar_type_SN == BLACK_HOLE & SNe_supernova_state == 2 & SNe_unbound == 1) | ...
                        (SNe_stellar_type_CP == BLACK_HOLE & SNe_supernova_state == 3 & SNe_unbound == 1));
number_three    = length(index_three);
data_three      = [SNe_mass_SN(index_three), SNe_component_speed_SN(index_three)];
 
h_three              = histogram(data_three(:,1),h_zero_bin_edges);
h_three_values       = h_three.Values;
h_three_bin_center   = h_three.BinEdges(1:end-1)+(h_three.BinWidth/2);

% 4. second born NS (when disrupted)
index_four      = find( (SNe_stellar_type_SN == NEUTRON_STAR & SNe_supernova_state == 2 & SNe_unbound == 1) | ...
                        (SNe_stellar_type_CP == NEUTRON_STAR & SNe_supernova_state == 3 & SNe_unbound == 1));
number_four     = length(index_four);
data_four       = [SNe_mass_SN(index_four), SNe_component_speed_SN(index_four)];

if debugFlag
    number_zero
    number_one
    number_two
    number_three
    number_four
end

if plotFlag
    % Plot
    clf
    figure()
    fs=18;
    lw=2.0;
    solar=char(9737);
    string1=['BH Mass [M_',solar,']'];  
    string2='Number of Systems';

    if metallicity==0.02
        Xlim = [2 20];
        text_x_lim = 10;
    else
        Xlim = [2 45];
        text_x_lim = 20;        
    end

    color0 = 105.*[1 1 1]./255;
    color1 = [0    0.4470    0.7410];
    color2 = [    0.8500    0.3250    0.0980];
    
    clf
    xlim(Xlim)
    xlabel(string1,'FontSize',fs)
    ylabel(string2,'FontSize',fs)
    ax=gca;
    ax.FontSize=fs;
    ax.FontName = 'Times New Roman';
    box on
    
    hold on
    plot(h_zero_bin_center,h_zero_values,'LineWidth',lw,'Color',color0)
    plot(h_one_bin_center,h_one_values,'LineWidth',lw,'Color',color1)
    plot(h_three_bin_center,h_three_values,'LineWidth',lw,'Color',color2)

    ylim=get(gca,'ylim');
    
    text(text_x_lim,ylim(2)-2500,strcat('Z=',num2str(metallicity)),'FontSize',fs,'FontName','Times New Roman')    
    
    legend( 'BH_{sin}',...
            'BH_1',...
            'BH_2',...
            'box','off',...
            'location','northeast',...
            'FontName','Times New Roman')
end

if savePlotFlag
        print(gcf,strcat('../plots/png/figure6_mass_distribution_Z=',num2str(metallicity),'.png'),'-dpng','-r300');
        saveas(gcf,strcat('../plots/pdf/figure6_mass_distribution_Z=',num2str(metallicity),'.pdf'))
end

end