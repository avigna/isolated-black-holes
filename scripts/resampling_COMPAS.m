function [  sys_realMassZAMS, ...
            isolated_real_total_mass_ZAMS,...
            isolated_kick_magnitude,...
            isolated_mass_SN,...
            isolated_stellar_type_SN] = resampling_COMPAS(filename,debugFlag,lowerLimit)

% MACROS
massEvolvedFactor = 4.0;    % Multiplicative factor to account for the non-simulated mass
                            % Obtained from integrating the full Kroupa IMF

% /BSE_System_Parameters
sys_massZAMS1           = h5read(filename,'/BSE_System_Parameters/Mass@ZAMS(1)');
sys_massZAMS2           = h5read(filename,'/BSE_System_Parameters/Mass@ZAMS(2)');
sys_SEED                = h5read(filename,'/BSE_System_Parameters/SEED');
sys_semi_major_axis     = h5read(filename,'/BSE_System_Parameters/Semi-Major_Axis');

% Calculate quantities
sys_totalMassZAMS       = sum(sys_massZAMS1)+sum(sys_massZAMS2);
sys_realMassZAMS        = massEvolvedFactor*sys_totalMassZAMS;

% Find primaries from wide binaries to use as single stars  
isolated_index                  = find(sys_semi_major_axis > lowerLimit);
isolated_SEEDs                  = sys_SEED(isolated_index);
isolated_total_mass_ZAMS        = sum(sys_massZAMS1(isolated_index))+sum(sys_massZAMS2(isolated_index));
isolated_real_total_mass_ZAMS   = massEvolvedFactor*isolated_total_mass_ZAMS;

% /BSE_Supernovae
SNe_applied_kick_magnitude  = h5read(filename,'/BSE_Supernovae/Applied_Kick_Magnitude(SN)');
SNe_mass_SN                 = h5read(filename,'/BSE_Supernovae/Mass(SN)');
SNe_SEED                    = h5read(filename,'/BSE_Supernovae/SEED');
SNe_stellar_type_SN         = h5read(filename,'/BSE_Supernovae/Stellar_Type(SN)');

% Look for the supernovae in the resampled population
[Lia,Locb] = ismember(SNe_SEED,isolated_SEEDs);
SNe_indexOfInterest = find(Lia==1);

isolated_kick_magnitude = SNe_applied_kick_magnitude(SNe_indexOfInterest);
isolated_mass_SN = SNe_mass_SN(SNe_indexOfInterest);
isolated_stellar_type_SN = SNe_stellar_type_SN(SNe_indexOfInterest);

if debugFlag
    display('Lower limit in sampling the separation [AU]:')
    lowerLimit
    display('Number of resampled systems:')
    length(isolated_index)
    display('Binary fraction:')
    length(sys_massZAMS1)/(length(sys_massZAMS1)+length(isolated_index))
end

end