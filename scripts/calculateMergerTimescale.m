function [T_Myr] = calculateMergerTimescale(M1,M2,a0,e0,capFlag)
% M1 is the mass of the primary in Msun
% M2 is the mass of the secondary in Msun
% a0 is the inital binary separation in AU
% e0 is the inital eccentricity
% The function returns T_Myr, the merger time, in Myr

% Based on Mandel 2021 fit
% https://arxiv.org/abs/2110.09254v1

Msunkg=1.98892e30;	%Msun in kg
c=299792458;		%speed of light, m/s
G=6.67428e-11;		%G in m^3/kg/s^2
AUtoRsun=215.032;   %AU to Rsun
RsolToM=6.957e+8;   %Rsun in meters
yearInSeconds = 3.154e+7;
AgeOfUniverseMyr = 13.772*1000;

Tc_Myr=(5/256).*(c^5.*(a0.*RsolToM*AUtoRsun).^4)./(G^3*(M1.*M2.*(M1+M2)*Msunkg^3))./yearInSeconds/10^6;
fe=(1+(0.27.*e0.^10)+(0.33.*e0.^20)+(0.2*e0.^1000)).*(1.0-(e0.^2)).^(7.0/2);

T_Myr = Tc_Myr.*fe;

if capFlag == 1
    T_Myr(find(T_Myr>AgeOfUniverseMyr))=AgeOfUniverseMyr;
end

end

