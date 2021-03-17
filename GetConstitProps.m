function [constitprops] = GetConstitProps
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
% Purpose: Specifies properties of composite constituent materials
% Input: None
% Output:
% - constitprops: Cell/struct containing constituent properties
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Preallocate
NconstitMax = 20;
constitprops = cell(1,NconstitMax);

%========================================================================
% -- Define all constituent materials below
%========================================================================

% -- Constituent Mat#1 - IM7 (Transversely Isotropic) 
Mat = 1;
constitprops{Mat}.name = "IM7";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 262.2E3; 
constitprops{Mat}.ET = 11.8E3;
constitprops{Mat}.vL = 0.17;
constitprops{Mat}.vT = 0.21;
constitprops{Mat}.GL = 18.9E3;
constitprops{Mat}.aL = -0.9E-06;
constitprops{Mat}.aT = 9.E-06;
constitprops{Mat}.allowables.XT = 4335.0;
constitprops{Mat}.allowables.XC = -2608.;
constitprops{Mat}.allowables.YT = 113.;
constitprops{Mat}.allowables.YC = -354.;
constitprops{Mat}.allowables.ZT = 113.;
constitprops{Mat}.allowables.ZC = -354.;
constitprops{Mat}.allowables.Q = 128.;
constitprops{Mat}.allowables.R = 138.;
constitprops{Mat}.allowables.S = 138.;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
GT = constitprops{Mat}.ET/(2*(1+constitprops{Mat}.vT));
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/GT;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#2 - Glass (Isotropic)
Mat = 2;
constitprops{Mat}.name = "Glass";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 73.E3; 
constitprops{Mat}.ET = 73.E3;
constitprops{Mat}.vL = 0.22;
constitprops{Mat}.vT = 0.22;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic  
constitprops{Mat}.aL = 5.E-06;
constitprops{Mat}.aT = 5.E-06;
constitprops{Mat}.allowables.XT = 2358.;
constitprops{Mat}.allowables.XC = -1653.;
constitprops{Mat}.allowables.YT = constitprops{Mat}.allowables.XT;
constitprops{Mat}.allowables.YC = constitprops{Mat}.allowables.XC;
constitprops{Mat}.allowables.ZT = constitprops{Mat}.allowables.XT;
constitprops{Mat}.allowables.ZC = constitprops{Mat}.allowables.XC;
constitprops{Mat}.allowables.Q = 1000.;
constitprops{Mat}.allowables.R = constitprops{Mat}.allowables.Q;
constitprops{Mat}.allowables.S = constitprops{Mat}.allowables.Q;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
GT = constitprops{Mat}.ET/(2*(1+constitprops{Mat}.vT));
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/GT;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#3 - 8552 (Isotropic)
Mat = 3;
constitprops{Mat}.name = "8552";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 4.67E3; 
constitprops{Mat}.ET = 4.67E3;
constitprops{Mat}.vL = 0.45;
constitprops{Mat}.vT = 0.45;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic  
constitprops{Mat}.aL = 42.E-06;
constitprops{Mat}.aT = 42.E-06;
AAA = 59.4; % -- Baseline
% AAA = 65; % -- MT
% AAA = 80.; % -- MOC
constitprops{Mat}.allowables.XT = AAA;
constitprops{Mat}.allowables.XC = -259.;
% constitprops{Mat}.allowables.XT = 1.E6; % -- High value to preclude axial failure
% constitprops{Mat}.allowables.XC = -1.E6;
constitprops{Mat}.allowables.YT = AAA;
constitprops{Mat}.allowables.YC = -259.;
constitprops{Mat}.allowables.ZT = AAA;
constitprops{Mat}.allowables.ZC = -259.;
constitprops{Mat}.allowables.Q = 112.0;
constitprops{Mat}.allowables.R = 112.0;
constitprops{Mat}.allowables.S = 112.0;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#4 - Epoxy (Isotropic) - Allowables based on MT
Mat = 4;
constitprops{Mat}.name = "Epoxy";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 3.45E3;  
constitprops{Mat}.ET = 3.45E3;   
constitprops{Mat}.vL = 0.35;
constitprops{Mat}.vT = 0.35;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic    
constitprops{Mat}.aL = 54.e-06; 
constitprops{Mat}.aT = 54.e-06; 
% -- Backed out based on MT - Tsai-Wu
AAA = 42.; 
BBB = -181.;
CCC = 81.;
constitprops{Mat}.allowables.XT = 1.E6;
constitprops{Mat}.allowables.XC = -1.E6;
constitprops{Mat}.allowables.YT = AAA;
constitprops{Mat}.allowables.YC = BBB;
constitprops{Mat}.allowables.ZT = AAA;
constitprops{Mat}.allowables.ZC = BBB;
constitprops{Mat}.allowables.Q = CCC;
constitprops{Mat}.allowables.R = CCC;
constitprops{Mat}.allowables.S = CCC;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#5 - Tsai&Hahn Glass (Isotropic)
Mat = 5;
constitprops{Mat}.name = "Tsai&Hahn Glass";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 113.4E3;  
constitprops{Mat}.ET = 113.4E3;   
constitprops{Mat}.vL = 0.22;
constitprops{Mat}.vT = 0.22;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic    
constitprops{Mat}.aL = 5.e-06; 
constitprops{Mat}.aT = 5.e-06; 

% -- Constituent Mat#6 - Tsai&Hahn Epoxy (Isotropic)
Mat = 6;
constitprops{Mat}.name = "Tsai&Hahn Epoxy";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 5.35E3;  
constitprops{Mat}.ET = 5.35E3;   
constitprops{Mat}.vL = 0.35;
constitprops{Mat}.vT = 0.35;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic    
constitprops{Mat}.aL = 54.e-06; 
constitprops{Mat}.aT = 54.e-06; 

% -- Constituent Mat#7 - Dean & Turner Carbon (Transversely Isotropic)
Mat = 7;
constitprops{Mat}.name = "Dean & Turner Carbon";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 232.0E3;  
constitprops{Mat}.ET = 15.0E3;   
constitprops{Mat}.vL = 0.279;
constitprops{Mat}.vT = 0.49;
constitprops{Mat}.GL = 24.E3;
constitprops{Mat}.aL = 0.e-06; % -- CTEs not given by Dean & Turner
constitprops{Mat}.aT = 0.e-06; 

% -- Constituent Mat#8 - Contour Glass (Isotropic)
Mat = 8;
constitprops{Mat}.name = "Contour Glass";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 69.E3;  
constitprops{Mat}.ET = 69.E3;   
constitprops{Mat}.vL = 0.2;
constitprops{Mat}.vT = 0.2;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic  
constitprops{Mat}.aL = 5.e-06; 
constitprops{Mat}.aT = 5.e-06; 

% -- Constituent Mat#9 - Contour Epoxy (Isotropic)
Mat = 9;
constitprops{Mat}.name = "Contour Epoxy";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 3.42E3;  
constitprops{Mat}.ET = 3.42E3;   
constitprops{Mat}.vL = 0.34;
constitprops{Mat}.vT = 0.34;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic   
constitprops{Mat}.aL = 54.e-06; 
constitprops{Mat}.aT = 54.e-06; 

% -- Constituent Mat#10 - SiC Fiber
Mat = 10;
constitprops{Mat}.name = "SiC Fiber";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 385.E3;  
constitprops{Mat}.ET = 385.E3;   
constitprops{Mat}.vL = 0.17;
constitprops{Mat}.vT = 0.17;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic    
constitprops{Mat}.aL = 3.2e-06; 
constitprops{Mat}.aT = 3.2e-06; 
constitprops{Mat}.allowables.XT = 1800.;  
constitprops{Mat}.allowables.XC = -1800.;  
constitprops{Mat}.allowables.YT = 1800.;
constitprops{Mat}.allowables.YC = -1800.;
constitprops{Mat}.allowables.ZT = 1800.;
constitprops{Mat}.allowables.ZC = -1800.;
constitprops{Mat}.allowables.Q = 700.;
constitprops{Mat}.allowables.R = 700.;
constitprops{Mat}.allowables.S = 700.;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
GT = constitprops{Mat}.ET/(2*(1+constitprops{Mat}.vT));
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/GT;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#11 - SiC Matrix
Mat = 11;
constitprops{Mat}.name = "SiC Matrix";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 327.E3;  
constitprops{Mat}.ET = 327.E3;   
constitprops{Mat}.vL = 0.22;
constitprops{Mat}.vT = 0.22;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL));   
constitprops{Mat}.aL = 3.1e-06; 
constitprops{Mat}.aT = 3.1e-06; 
constitprops{Mat}.allowables.XT = 600.;  
constitprops{Mat}.allowables.XC = -600.;  
constitprops{Mat}.allowables.YT = 600.;
constitprops{Mat}.allowables.YC = -600.;
constitprops{Mat}.allowables.ZT = 600.;
constitprops{Mat}.allowables.ZC = -600.;
constitprops{Mat}.allowables.Q = 350.;
constitprops{Mat}.allowables.R = 350.;
constitprops{Mat}.allowables.S = 350.;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
GT = constitprops{Mat}.ET/(2*(1+constitprops{Mat}.vT));
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/GT;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#12 - BN Coating
Mat = 12;
constitprops{Mat}.name = "BN Coating";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 10.E3;  
constitprops{Mat}.ET = 10.E3;   
constitprops{Mat}.vL = 0.23;
constitprops{Mat}.vT = 0.23;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic    
constitprops{Mat}.aL = 4.e-06; 
constitprops{Mat}.aT = 4.e-06; 
X = 70.;
Q = X*45/70;
constitprops{Mat}.allowables.XT = X;  
constitprops{Mat}.allowables.XC = -X;  
constitprops{Mat}.allowables.YT = X;
constitprops{Mat}.allowables.YC = -X;
constitprops{Mat}.allowables.ZT = X;
constitprops{Mat}.allowables.ZC = -X;
constitprops{Mat}.allowables.Q = Q;
constitprops{Mat}.allowables.R = Q;
constitprops{Mat}.allowables.S = Q;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
GT = constitprops{Mat}.ET/(2*(1+constitprops{Mat}.vT));
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/GT;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#13 - Sun & Vaidya Carbon
Mat = 13;
constitprops{Mat}.name = "Sun & Vaidya Carbon";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 235.E3;  
constitprops{Mat}.ET = 14.E3;   
constitprops{Mat}.vL = 0.2;
constitprops{Mat}.vT = 0.25;
constitprops{Mat}.GL = 28.E3;   
constitprops{Mat}.aL = -0.9e-06; 
constitprops{Mat}.aT = 9.e-06; 

% -- Constituent Mat#14 - Sun & Vaidya Epoxy (Isotropic)
Mat = 14;
constitprops{Mat}.name = "Sun & Vaidya Epoxy";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 4.8E3;  
constitprops{Mat}.ET = 4.8E3;   
constitprops{Mat}.vL = 0.34;
constitprops{Mat}.vT = 0.34;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic   
constitprops{Mat}.aL = 54.e-06; 
constitprops{Mat}.aT = 54.e-06; 

% -- Constituent Mat#15 - Epoxy (Isotropic) - backed out from HF - Max Stress - 46x80 hex
Mat = 15;
constitprops{Mat}.name = "Epoxy - HF - Max Stress - 46x80 hex";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 3.45E3;  
constitprops{Mat}.ET = 3.45E3;   
constitprops{Mat}.vL = 0.35;
constitprops{Mat}.vT = 0.35;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL));   
constitprops{Mat}.aL = 54.e-06; 
constitprops{Mat}.aT = 54.e-06; 
%----------- Chapter 7 backed out strengths ------------
TTT = 86.;
CCC = -425;
SSS = 184;
constitprops{Mat}.allowables.XT = TTT;
constitprops{Mat}.allowables.XC = CCC;
constitprops{Mat}.allowables.YT = TTT;
constitprops{Mat}.allowables.YC = CCC;
constitprops{Mat}.allowables.ZT = TTT;
constitprops{Mat}.allowables.ZC = CCC;
constitprops{Mat}.allowables.Q = SSS;
constitprops{Mat}.allowables.R = SSS;
constitprops{Mat}.allowables.S = SSS;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#16 - damaged Epoxy (Isotropic) from Mat#15
Mat = 16;
constitprops{Mat}.name = "Damaged Epoxy";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 3.45E3/10000.;  
constitprops{Mat}.ET = 3.45E3/10000;   
constitprops{Mat}.vL = 0.35;
constitprops{Mat}.vT = 0.35;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL));   
constitprops{Mat}.aL = 54.e-06; 
constitprops{Mat}.aT = 54.e-06; 
%----------- Already damaged, so strengths set very high ------------
TTT = 86.E6;
CCC = -425E6;
SSS = 184E6;
constitprops{Mat}.allowables.XT = TTT;
constitprops{Mat}.allowables.XC = CCC;
constitprops{Mat}.allowables.YT = TTT;
constitprops{Mat}.allowables.YC = CCC;
constitprops{Mat}.allowables.ZT = TTT;
constitprops{Mat}.allowables.ZC = CCC;
constitprops{Mat}.allowables.Q = SSS;
constitprops{Mat}.allowables.R = SSS;
constitprops{Mat}.allowables.S = SSS;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

% -- Constituent Mat#17 - BN Coating with Failure OFF
Mat = 17;
constitprops{Mat}.name = "BN Coating";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 10.E3;  
constitprops{Mat}.ET = 10.E3;   
constitprops{Mat}.vL = 0.23;
constitprops{Mat}.vT = 0.23;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic    
constitprops{Mat}.aL = 4.e-06; 
constitprops{Mat}.aT = 4.e-06; 

%Constituent Mat#18 - Glass (Isotropic) - backed out from HF - Max Stress - 46x80 hex
Mat = 18;
constitprops{Mat}.name = "Glass - HF - Max Stress - 46x80 hex";
constitprops{Mat}.constID = Mat;
constitprops{Mat}.EL = 73.E3; 
constitprops{Mat}.ET = 73.E3;
constitprops{Mat}.vL = 0.22;
constitprops{Mat}.vT = 0.22;
constitprops{Mat}.GL = constitprops{Mat}.EL/(2*(1+constitprops{Mat}.vL)); % -- isotropic  
constitprops{Mat}.aL = 5.E-06;
constitprops{Mat}.aT = 5.E-06;
% -- Chapter 7 backed out strengths
constitprops{Mat}.allowables.XT = 2445.;
constitprops{Mat}.allowables.XC = -1653.;
constitprops{Mat}.allowables.YT = constitprops{Mat}.allowables.XT;
constitprops{Mat}.allowables.YC = constitprops{Mat}.allowables.XC;
constitprops{Mat}.allowables.ZT = constitprops{Mat}.allowables.XT;
constitprops{Mat}.allowables.ZC = constitprops{Mat}.allowables.XC;
constitprops{Mat}.allowables.Q = 1000.;
constitprops{Mat}.allowables.R = constitprops{Mat}.allowables.Q;
constitprops{Mat}.allowables.S = constitprops{Mat}.allowables.Q;
constitprops{Mat}.allowables.XeT = constitprops{Mat}.allowables.XT/constitprops{Mat}.EL;
constitprops{Mat}.allowables.XeC = constitprops{Mat}.allowables.XC/constitprops{Mat}.EL;
constitprops{Mat}.allowables.YeT = constitprops{Mat}.allowables.YT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.YeC = constitprops{Mat}.allowables.YC/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeT = constitprops{Mat}.allowables.ZT/constitprops{Mat}.ET;
constitprops{Mat}.allowables.ZeC = constitprops{Mat}.allowables.ZC/constitprops{Mat}.ET;
GT = constitprops{Mat}.ET/(2*(1+constitprops{Mat}.vT));
constitprops{Mat}.allowables.Qe = constitprops{Mat}.allowables.Q/GT;
constitprops{Mat}.allowables.Re = constitprops{Mat}.allowables.R/constitprops{Mat}.GL;
constitprops{Mat}.allowables.Se = constitprops{Mat}.allowables.S/constitprops{Mat}.GL;

%========================================================================
% -- End of constituent material definitions
%========================================================================

% -- Ensure that if a material has allowables, it has all of them
for Mat = 1: length(constitprops)
    if isfield(constitprops{Mat}, 'allowables')
        Names = {"XT", "XC", "YT", "YC", "ZT", "ZC", "Q", "R", "S", ...
                 "XeT","XeC","YeT","YeC","ZeT","ZeC","Qe","Re","Se"};
        for i = 1: length(Names)
            [constitprops{Mat}.allowables] = ...
                SetUndefinedAllowable(constitprops{Mat}.allowables, Names{i});
        end
    end
end

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [allowables] = SetUndefinedAllowable(allowables, Name)
    if ~isfield(allowables, char(Name))
        if strfind(Name,"C") > 1
            allowables.(Name) = -1.E99;
        else
            allowables.(Name) = 1.E99;
        end
    end
end


