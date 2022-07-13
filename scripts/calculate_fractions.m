function calculate_fractions(filename, label, debugFlag, savaDataFlag, lowerLimit)

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
    number_two
    number_three
    number_four
end

% Gather and save data
% Single
single_fraction_mass_2_5      = calculate_single_fraction(2.5,data_zero,data_one,data_three);
single_fraction_mass_5        = calculate_single_fraction(5,data_zero,data_one,data_three);
single_fraction_mass_10       = calculate_single_fraction(10,data_zero,data_one,data_three);
single_fraction_mass_15       = calculate_single_fraction(15,data_zero,data_one,data_three);
single_fraction_mass_20       = calculate_single_fraction(20,data_zero,data_one,data_three);
single_fraction_mass_100      = calculate_single_fraction(100,data_zero,data_one,data_three);

% Binary
binary_fraction_mass_2_5      = calculate_binary_fraction(2.5,data_zero,data_one,data_three);
binary_fraction_mass_5        = calculate_binary_fraction(5,data_zero,data_one,data_three);
binary_fraction_mass_10       = calculate_binary_fraction(10,data_zero,data_one,data_three);
binary_fraction_mass_15       = calculate_binary_fraction(15,data_zero,data_one,data_three);
binary_fraction_mass_20       = calculate_binary_fraction(20,data_zero,data_one,data_three);
binary_fraction_mass_100      = calculate_binary_fraction(100,data_zero,data_one,data_three);

% Single fraction (larger than given mass)
fraction_abs_BH_sin_mass_2_5    = calculate_single_absolute_fraction(2.5,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_4      = calculate_single_absolute_fraction(4,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_5      = calculate_single_absolute_fraction(5,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_6      = calculate_single_absolute_fraction(6,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_8      = calculate_single_absolute_fraction(8,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_10     = calculate_single_absolute_fraction(10,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_12     = calculate_single_absolute_fraction(12,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_14     = calculate_single_absolute_fraction(14,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_15     = calculate_single_absolute_fraction(15,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_16     = calculate_single_absolute_fraction(16,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_18     = calculate_single_absolute_fraction(18,data_zero,data_one,data_three);
fraction_abs_BH_sin_mass_20     = calculate_single_absolute_fraction(20,data_zero,data_one,data_three);


if debugFlag
    single_fraction_mass_5
    single_fraction_mass_10
    single_fraction_mass_15
    single_fraction_mass_20
    single_fraction_mass_100

    binary_fraction_mass_5
    binary_fraction_mass_10
    binary_fraction_mass_15
    binary_fraction_mass_20
    binary_fraction_mass_100

    fraction_abs_BH_sin_mass_2_5
    fraction_abs_BH_sin_mass_4
    fraction_abs_BH_sin_mass_5
    fraction_abs_BH_sin_mass_6
    fraction_abs_BH_sin_mass_8
    fraction_abs_BH_sin_mass_10
    fraction_abs_BH_sin_mass_12
    fraction_abs_BH_sin_mass_14
    fraction_abs_BH_sin_mass_15
    fraction_abs_BH_sin_mass_16
    fraction_abs_BH_sin_mass_18
    fraction_abs_BH_sin_mass_20
end


if savaDataFlag
    save(   strcat('../data/fraction_Z_',label,'.mat'),...
            'single_fraction_mass_2_5',...
            'single_fraction_mass_5',...
            'single_fraction_mass_10',...
            'single_fraction_mass_15',...
            'single_fraction_mass_20',...
            'single_fraction_mass_100',...
            'binary_fraction_mass_2_5',...
            'binary_fraction_mass_5',...
            'binary_fraction_mass_10',...
            'binary_fraction_mass_15',...
            'binary_fraction_mass_20',...
            'binary_fraction_mass_100',...
            'fraction_abs_BH_sin_mass_2_5',...
            'fraction_abs_BH_sin_mass_4',...
            'fraction_abs_BH_sin_mass_5',...
            'fraction_abs_BH_sin_mass_6',...
            'fraction_abs_BH_sin_mass_8',...
            'fraction_abs_BH_sin_mass_10',...
            'fraction_abs_BH_sin_mass_12',...
            'fraction_abs_BH_sin_mass_14',...
            'fraction_abs_BH_sin_mass_15',...
            'fraction_abs_BH_sin_mass_16',...
            'fraction_abs_BH_sin_mass_18',...
            'fraction_abs_BH_sin_mass_20')
end

end

function [fraction_BH_bin] =calculate_binary_fraction(mass_limit,data_zero,data_one,data_three)
    N_BH_sin    = nnz(find((data_zero(:,1) <= mass_limit)));
    N_BH_1      = nnz(find((data_one(:,1) <= mass_limit)));
    N_BH_2      = nnz(find((data_three(:,1) <= mass_limit)));
    N_BH_bin    = N_BH_1 + N_BH_2;
    fraction_BH_bin = N_BH_bin/(N_BH_bin+N_BH_sin);
end

function [fraction_BH_sin] =calculate_single_fraction(mass_limit,data_zero,data_one,data_three)
    N_BH_sin    = nnz(find((data_zero(:,1) <= mass_limit)));
    N_BH_1      = nnz(find((data_one(:,1) <= mass_limit)));
    N_BH_2      = nnz(find((data_three(:,1) <= mass_limit)));
    N_BH_bin    = N_BH_1 + N_BH_2;
    fraction_BH_sin = N_BH_sin/(N_BH_bin+N_BH_sin);
end

function [fraction_abs_BH_sin] =calculate_single_absolute_fraction(mass_limit,data_zero,data_one,data_three)
    N_BH_sin    = nnz(find((data_zero(:,1) >= mass_limit)));
    N_BH_1      = nnz(find((data_one(:,1) >= mass_limit)));
    N_BH_2      = nnz(find((data_three(:,1) >= mass_limit)));
    N_BH_bin    = N_BH_1 + N_BH_2;
    fraction_abs_BH_sin = N_BH_sin/(N_BH_bin+N_BH_sin);
end