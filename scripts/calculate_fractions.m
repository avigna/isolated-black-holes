function calculate_fractions(filename, label, debugFlag, saveFlag)

% MACROS
% Stellar types as used in COMPAS, originally defined by Hurley+2000
NEUTRON_STAR            = 13;
BLACK_HOLE              = 14;

% Get isolated systems
lowerLimit = 80;
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

% % Calculate quantities
% SNe_total_mass      = SNe_mass_CP + SNe_mass_SN;
% SNe_periapsis       = SNe_SemiMajorAxis.*(1-SNe_eccentricity);

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
 

% 4. second born NS (when disrupted)
index_four      = find( (SNe_stellar_type_SN == NEUTRON_STAR & SNe_supernova_state == 2 & SNe_unbound == 1) | ...
                        (SNe_stellar_type_CP == NEUTRON_STAR & SNe_supernova_state == 3 & SNe_unbound == 1));
number_four     = length(index_four);
data_four       = [SNe_mass_SN(index_four), SNe_component_speed_SN(index_four)];

if debugFlag
    number_zero
    number_one
    number_three
end

% Gather and save data
fraction_mass_5_speed_30        = calculate_binary_fraction(5,30,data_zero,data_one,data_three);
fraction_mass_10_speed_30       = calculate_binary_fraction(10,30,data_zero,data_one,data_three);
fraction_mass_20_speed_30       = calculate_binary_fraction(20,30,data_zero,data_one,data_three);
fraction_mass_100_speed_30      = calculate_binary_fraction(100,30,data_zero,data_one,data_three);
fraction_mass_5_speed_50        = calculate_binary_fraction(5,50,data_zero,data_one,data_three);
fraction_mass_10_speed_50       = calculate_binary_fraction(10,50,data_zero,data_one,data_three);
fraction_mass_20_speed_50       = calculate_binary_fraction(20,50,data_zero,data_one,data_three);
fraction_mass_100_speed_50      = calculate_binary_fraction(100,50,data_zero,data_one,data_three);
fraction_mass_5_speed_100       = calculate_binary_fraction(5,100,data_zero,data_one,data_three);
fraction_mass_10_speed_100      = calculate_binary_fraction(10,100,data_zero,data_one,data_three);
fraction_mass_20_speed_100      = calculate_binary_fraction(20,100,data_zero,data_one,data_three);
fraction_mass_100_speed_100      = calculate_binary_fraction(100,100,data_zero,data_one,data_three);

fraction_mass_10_speed_50_greater_than =calculate_binary_fraction_greater_than(10,50,data_zero,data_one,data_three);
fraction_mass_20_speed_50_greater_than =calculate_binary_fraction_greater_than(20,50,data_zero,data_one,data_three);

if debugFlag
    fraction_mass_5_speed_30
    fraction_mass_10_speed_30
    fraction_mass_20_speed_30
    fraction_mass_100_speed_30
    fraction_mass_5_speed_50
    fraction_mass_10_speed_50
    fraction_mass_20_speed_50
    fraction_mass_100_speed_50
    fraction_mass_5_speed_100
    fraction_mass_10_speed_100
    fraction_mass_20_speed_100
    fraction_mass_100_speed_100
    fraction_mass_10_speed_50_greater_than
    fraction_mass_20_speed_50_greater_than
end


if saveFlag
    save(   strcat('../data/fraction_Z_',label,'.mat'),...
            'fraction_mass_5_speed_30',...
            'fraction_mass_10_speed_30',...
            'fraction_mass_20_speed_30',...
            'fraction_mass_100_speed_30',...
            'fraction_mass_5_speed_50',...
            'fraction_mass_10_speed_50',...
            'fraction_mass_20_speed_50',...
            'fraction_mass_100_speed_50',...            
            'fraction_mass_5_speed_100',...
            'fraction_mass_10_speed_100',...
            'fraction_mass_20_speed_100',...
            'fraction_mass_100_speed_100')
end

end

function [fraction_BH_bin] =calculate_binary_fraction(mass_limit,velocity_limit,data_zero,data_one,data_three)
    N_BH_sin    = nnz(find((data_zero(:,1) <= mass_limit) & (data_zero(:,2) <= velocity_limit)));
    N_BH_1      = nnz(find((data_one(:,1) <= mass_limit) & (data_one(:,2) <= velocity_limit)));
    N_BH_2      = nnz(find((data_three(:,1) <= mass_limit) & (data_three(:,2) <= velocity_limit)));
    N_BH_bin    = N_BH_1 + N_BH_2;
    fraction_BH_bin = N_BH_bin/(N_BH_bin+N_BH_sin);
end

function [fraction_BH_bin_greaterthan] =calculate_binary_fraction_greater_than(mass_limit,velocity_limit,data_zero,data_one,data_three)
    N_BH_sin    = nnz(find((data_zero(:,1) > mass_limit) & (data_zero(:,2) <= velocity_limit)));
    N_BH_1      = nnz(find((data_one(:,1) > mass_limit) & (data_one(:,2) <= velocity_limit)));
    N_BH_2      = nnz(find((data_three(:,1) > mass_limit) & (data_three(:,2) <= velocity_limit)));
    N_BH_bin    = N_BH_1 + N_BH_2;
    fraction_BH_bin_greaterthan = N_BH_bin/(N_BH_bin+N_BH_sin);
end