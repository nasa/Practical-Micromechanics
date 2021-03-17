function [MoS] = FCriteria(stress, strain, stressth, strainth, allowables, ...
                           CriteriaOn, micro, MoS)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Copyright 2020 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration. No copyright is claimed in the 
% United States under Title 17, U.S. Code. All Other Rights Reserved. BY DOWNLOADING 
% OR USING THIS SOFTWARE, YOU ACKNOWLEDGE THAT YOU HAVE READ THE NASA OPEN SOURCE 
% AGREEMENT V1.3, THAT YOU UNDERSTAND IT, AND THAT YOU AGREE TO BE BOUND BY ITS 
% TERMS. IF YOU DO NOT AGREE TO THE TERMS AND CONDITIONS OF THIS AGREEMENT, DO NOT 
% USE OR DOWNLOAD THE SOFTWARE. THIS SOFTWARE IS PROVIDED AS IS WITHOUT ANY WARRANTY 
% OF ANY KIND. RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST, AND INDEMNIFIES 
% AND HOLDS HARMLESS, THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND 
% SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT. This code was prepared by Drs. 
% B.A. Bednarcyk and S.M. Arnold to complement the book “Practical Micromechanics of 
% Composite Materials” during the course of their government work.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Given a stress and strain state, calculate the MoS by calling functions 
%          (contained herein) for each failure criterion 
% Input:
% - stress: Stresses (from mechanical problem only if ThermoMech) 
% - strain: Strains (from mechanical problem only if ThermoMech) 
% - stressth: Stresses from thermal problem only if ThermoMech, 0 otherwise
% - strainth: Strains from thermal problem only if ThermoMech, 0 otherwise
% - allowables: Struct containing material allowables
% - CriteriaOn: Flags indicating which failure criteria are turned on
% - micro: Flag indicating if problem is micromechanics-based
% - MoS: Struct containing MoS results and info
% Output:
% - MoS: Updated struct containing MoS results and info
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Determine if problem is combined thermomechanical
thermomech = false;
for i = 1:6
    if stressth(i) ~= 0
        thermomech = true;
    end
end

% -- Check if thermal preload causes failure
if thermomech
    zero6 = zeros(6,1);
    
    % -- Calculate margins from thermal preload only
    [MoSTh.Crit{1}] = MaxStress(stressth, zero6, allowables, micro);
    [MoSTh.Crit{2}] = MaxStrain(strainth, zero6, allowables, micro);
    [MoSTh.Crit{3}] = TsaiHill(stressth, zero6, allowables, micro);
    [MoSTh.Crit{4}] = TsaiWu(stressth, zero6, allowables, micro);
    
    for j = 1:4
       if j < 3
           num = 6;
       else
           num = 1;
       end
       for i = 1:num
           if MoSTh.Crit{j}.MoS(i) <= 0
                error(['** Negative Margin Detected Due to Thermal Preload ', ...
                       ' To see results, run with thermal load only **'])
           end
       end
    end

end

% -- Standard margin calculations, call each criterion
[MoS.Crit{1}] = MaxStress(stress, stressth, allowables, micro);
[MoS.Crit{2}] = MaxStrain(strain, strainth, allowables, micro);
[MoS.Crit{3}] = TsaiHill(stress, stressth, allowables, micro);
[MoS.Crit{4}] = TsaiWu(stress, stressth, allowables, micro);

% -- Set MoS high if criterion is turned off
for i = 1:4
    if CriteriaOn(i) == 0
        MoS.Crit{i}.MoS = [99999, 99999, 99999, 99999, 99999, 99999];
        MoS.Crit{i}.MinMoS = 99999;
    end
end

% -- Set and store controlling MoS information
MoS.Controlling = "none";
MoS.MinMoS = 99999;
MoS.ControllingNum = 0;
for i = 1:4
    if MoS.Crit{i}.MinMoS < MoS.MinMoS
        MoS.Controlling = MoS.Crit{i}.name;
        MoS.ControllingNum = i;
        MoS.MinMoS = MoS.Crit{i}.MinMoS;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Crit] = MaxStress(stress, stressth, allowables, micro)

% -- Max stress failure criterion

Crit.name = "MaxStress";

% -- Calculate stress ratios (R_i)
[Rmech, Rth, ~] = GetR(stress, stressth, allowables, micro);

% -- Calculate MoS and determine controlling component
Crit.Controlling = 0;
Crit.MinMoS = 99999;
for i = 1:6 % -- Loop through 6 stress component-based MoS
    if abs(Rmech(i)) < 0.0000001 % -- Set MoS high if R is small
        Crit.MoS(i) = 99999;
    else
        Crit.MoS(i) = (1 - Rth(i))/Rmech(i) - 1; % -- Eqs. 4.45 & 4.46 
    end
    if Crit.MoS(i) < Crit.MinMoS % -- Store controlling MoS
        Crit.MinMoS = Crit.MoS(i);
        Crit.Controlling = i;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Crit] = MaxStrain(strain, strainth, allowables, micro)

% -- Max strain failure criterion

% -- Note:
%    - strain is the strain from the mechanical load case
%    - strainth is the mechanical strain from the thermal load case

Crit.name = "MaxStrain";

strainallowables.XT = allowables.XeT;
strainallowables.XC = allowables.XeC;
strainallowables.YT = allowables.YeT;
strainallowables.YC = allowables.YeC;
strainallowables.S = allowables.Se;
if micro
    strainallowables.ZT = allowables.ZeT;
    strainallowables.ZC = allowables.ZeC;
    strainallowables.Q = allowables.Qe;
    strainallowables.R = allowables.Re;
end
    
% -- Calculate strain ratios (R_i)
[Rmech, Rth, ~] = GetR(strain, strainth, strainallowables, micro);

% -- Calculate MoS and determine controlling component
Crit.Controlling = 0;
Crit.MinMoS = 99999;
for i = 1:6 % -- Loop through 6 strain component-based MoS
    if abs(Rmech(i)) < 0.0000001 % -- Set MoS high if R is small
        Crit.MoS(i) = 99999;
    else
        Crit.MoS(i) = (1-Rth(i))/Rmech(i) - 1; % -- Eqs. 4.45 - 4.47
    end
    if Crit.MoS(i) < Crit.MinMoS % -- Store controlling MoS
        Crit.MinMoS = Crit.MoS(i);
        Crit.Controlling = i;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Crit] = TsaiHill(stress, stressth, allowables, micro)

% -- Tsai-Hill failure criterion

% -- Calculate stress ratios (R_i) and current allowables (based on T vs C)
[Rmech, Rth, CurrentAllowables] = GetR(stress, stressth, allowables, micro);

X = CurrentAllowables.X;
Y = CurrentAllowables.Y;
Z = CurrentAllowables.Z;

% -- Determine D and K coefficients (Eq. 4.50)
D = zeros(3,1);
K = eye(6,6);

K(2,3) = -Y*Z*(-1/X^2 + 1/Y^2 + 1/Z^2);
K(1,3) = -X*Z*(1/X^2 - 1/Y^2 + 1/Z^2);
K(1,2) = -X*Y*(1/X^2 + 1/Y^2 - 1/Z^2);

% -- Calculate MoS
[Crit] = TsaiHillWuMoS(Rmech, Rth, D, K);

Crit.name = "Tsai-Hill";

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Crit] = TsaiWu(stress, stressth, allowables, micro)

% -- Tsai-Wu failure criterion

% -- Calculate stress ratios (R_i) and current allowables (based on T vs C)
[Rmech, Rth, CurrentAllowables] = GetR(stress, stressth, allowables, micro);

X = CurrentAllowables.X;
Y = CurrentAllowables.Y;
Z = CurrentAllowables.Z;

XT = allowables.XT;
XC = allowables.XC;
YT = allowables.YT;
YC = allowables.YC;
S = allowables.S;

if micro % -- Micromechanics-based
    Q = allowables.Q;
    R = allowables.R;
    ZT = allowables.ZT;
    ZC = allowables.ZC;
else % -- Not micromechanics-based, Q,R,Z not used, but can't be 0
    Q = 1;
    R = 1;
    ZT = 1;
    ZC = 1;
end

% -- Eq 4.8
F1 = 1/XT + 1/XC;
F2 = 1/YT + 1/YC;
F3 = 1/ZT + 1/ZC;
F11 = -1/(XT*XC);
F22 = -1/(YT*YC);
F33 = -1/(ZT*ZC);
F44 = 1/Q^2;
F55 = 1/R^2;
F66 = 1/S^2;

% -- Tsai-Hahn (1980) interaction coefficients (Eq. 4.9)
F12 = -1/(2*sqrt(XT*XC*YT*YC));
F13 = -1/(2*sqrt(XT*XC*ZT*ZC));
F23 = -1/(2*sqrt(YT*YC*ZT*ZC));

% -- Determine D and K coefficients (Eq. 4.51)
D = zeros(3,1);
K = zeros(6,6);
D(1) = F1*X;
D(2) = F2*Y;
D(3) = F3*Z;
K(1,1) = F11*X^2;
K(2,2) = F22*Y^2;
K(3,3) = F33*Z^2;
K(4,4) = F44*Q^2;
K(5,5) = F55*R^2;
K(6,6) = F66*S^2;
K(2,3) = 2*F23*Y*Z;
K(1,3) = 2*F13*X*Z;
K(1,2) = 2*F12*X*Y;

% -- Calculate MoS
[Crit] = TsaiHillWuMoS(Rmech, Rth, D, K);

Crit.name = "Tsai-Wu";

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, Rth, CurrentAllowables] = GetR(Smech, Sth, allowables, micro)

% -- Calculate the stress ratios (R_i) and the current allowables (based on
%    tension vs. compression)

Rth = zeros(1,6);

if Smech(1) == 0
    R(1) = 0.00000001; % -- Rmech can't be zero for max stress & max strain
    Rth(1) = Sth(1)/allowables.XT;
    CurrentAllowables.X = allowables.XT;
elseif Smech(1) > 0
    R(1) = Smech(1)/allowables.XT;
    Rth(1) = Sth(1)/allowables.XT;
    CurrentAllowables.X = allowables.XT;
else
    R(1) = Smech(1)/allowables.XC;
    Rth(1) = Sth(1)/allowables.XC;
    CurrentAllowables.X = allowables.XC;
end

if Smech(2) == 0
    R(2) = 0.00000001;
    Rth(2) = Sth(2)/allowables.YT;
    CurrentAllowables.Y = allowables.YT;
elseif Smech(2) > 0
    R(2) = Smech(2)/allowables.YT;
    Rth(2) = Sth(2)/allowables.YT;
    CurrentAllowables.Y = allowables.YT;
else
    R(2) = Smech(2)/allowables.YC;
    Rth(2) = Sth(2)/allowables.YC;
    CurrentAllowables.Y = allowables.YC;
end

if Smech(6) == 0
    R(6) = 0.00000001;
else
    R(6) = abs(Smech(6))/allowables.S;
end 
Rth(6) = abs(Sth(6))/allowables.S;
CurrentAllowables.S = allowables.S;


if micro % -- Micromechanics-based (fully 3D)
    if Smech(3) == 0
        R(3) = 0.00000001;
        Rth(3) = Sth(3)/allowables.ZT;
        CurrentAllowables.Z = allowables.ZT;
    elseif Smech(3) > 0
        R(3) = Smech(3)/allowables.ZT;
        Rth(3) = Sth(3)/allowables.ZT;
        CurrentAllowables.Z = allowables.ZT;
    else
        R(3) = Smech(3)/allowables.ZC;
        Rth(3) = Sth(3)/allowables.ZC;
        CurrentAllowables.Z = allowables.ZC;
    end 
    
    if Smech(4) == 0
        R(4) = 0.00000001;
    else
        R(4) = abs(Smech(4))/allowables.Q;
    end
    Rth(4) = abs(Sth(4))/allowables.Q;
    CurrentAllowables.Q = allowables.Q;
   
    if Smech(5) == 0
        R(5) = 0.00000001;
    else
        R(5) = abs(Smech(5))/allowables.R;
    end   
    Rth(5) = abs(Sth(5))/allowables.R;
    CurrentAllowables.R = allowables.R;
    
else % -- Not micromechanics-based Q and R not used, Z = Y for Tsai-Hill
    CurrentAllowables.Q = 1;
    CurrentAllowables.R = 1;
    CurrentAllowables.Z = CurrentAllowables.Y;
        
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Crit] = TsaiHillWuMoS(Rmech, Rth, D, K)

% -- Calculate actual MoS for Tsai-Hill and Tsai-Wu

% -- Eq. 4.53
P1mech = Rmech(1:3)*D;
P1th = Rth(1:3)*D;

% -- Eq. 4.55
P2mech = 0;
P2mechth = 0;
P2th = 0;
for i = 1:6
    for j = 1:6
        P2mech = P2mech + K(i,j)*Rmech(i)*Rmech(j);
        P2mechth = P2mechth + K(i,j)*(Rmech(i)*Rth(j) + Rth(i)*Rmech(j));
        P2th = P2th + K(i,j)*Rth(i)*Rth(j);
    end
end

% -- Eq. 4.57
A = P2mech;
B = P2mechth + P1mech;
C = P2th + P1th - 1;

% -- Check for sqrt term in 4.58 less than zero (should only be numerics)
if B^2 - 4*A*C < 0
    den = 0;
else
    den = B + sqrt(B^2 - 4*A*C); 
end
 
if den == 0
    Crit.MoS = 99999;
else
    Crit.MoS = -2*C/den - 1; % -- Eq. 4.58
end

if Crit.MoS > 99999
    Crit.MoS = 99999;
end

% -- Set controlling to 1 - there is only one MoS for Tsai-Hill and Tsai-Wu
Crit.Controlling = 1;
Crit.MinMoS = Crit.MoS;

end