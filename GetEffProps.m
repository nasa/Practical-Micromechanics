function [props] = GetEffProps(constitprops, props)
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
%
% Purpose: Specifies material number, name, micromechanics theory, constituent volume
%          fractions, and constituent materials for micromechanics-based composite 
%          materials, and runs the micromechanics theory to obtain effective props 
%          and concentration tensors.
%          ** Any Mat not used in a Problem will not be evaluated **
% Input:
% - constitprops: Cell/struct containing constituent properties
% - props: Cell/struct containing effective composite properties
% Output:
% - props: Updated cell/struct containing effective composite properties
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('*** Evaluating Micromechanics-Based Materials ***');

%========================================================================
% -- Define all composite materials below
% -- *** Begin numbering at 100 ***
%========================================================================

Mat = 100;
name = "Voigt IM7-8552 Vf = 0.55";
Theory = 'Voigt';
Vf = 0.55;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 101;
name = "Reuss IM7-8552 Vf = 0.55";
Theory = 'Reuss';
Vf = 0.55;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 102;
name = "MT IM7-8552 Vf = 0.55";
Theory = 'MT';
Vf = 0.55;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 103;
name = "MOC IM7-8552 Vf = 0.55";
Theory = 'MOC';
Vf = 0.55;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 104;
name = "MOCu IM7-8552 Vf = 0.55";
Theory = 'MOCu';
Vf = 0.55;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 105;
name = "All theories IM7-8552 Vf = x";
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
% -- Call function (contained in this file) to get eff props across all Vf
RunOverVf(name, Constits, props, Mat);
props{Mat}.Quit = true; % -- Tells MicroAnalysis to quit after getting eff props

Mat = 106;
name = "All theories glass-epoxy Vf = x";
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
% -- Call function (contained in this file) to get eff props across all Vf
RunOverVf(name, Constits, props, Mat);
props{Mat}.Quit = true; % -- Tells MicroAnalysis to quit after getting eff props


Mat = 130;
name = "2x2-7x7-and-26x26 Glass-Epoxy Tsai-Hahn GMC - Fig 5.4";
Constits.Fiber = constitprops{5};
Constits.Matrix = constitprops{6};
% -- Call function (contained in this file) to get eff props across all Vf
RunOverVf(name, Constits, props, Mat);
props{Mat}.Quit = true; % -- Tells MicroAnalysis to quit after getting eff props

Mat = 131;
name = "2x2-7x7-and-26x26 Carbon-Epoxy Dean & Turner GMC - Fig 5.5";
Constits.Fiber = constitprops{7};
Constits.Matrix = constitprops{6};
% -- Call function (contained in this file) to get eff props across all Vf
RunOverVf(name, Constits, props, Mat);
props{Mat}.Quit = true; % -- Tells MicroAnalysis to quit after getting eff props

Mat = 132;
Theory = 'GMC';
name = "26x26 Glass-Epoxy Contour GMC Fig 5.10";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);
 
Mat = 133;
Theory = 'GMC';
name = "2x2 SiC-SiC GMC Table 5.3";
Vf = 0.28;
Vi  = 0;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
RUCid = 2;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);
 
Mat = 134;
Theory = 'GMC';
name = "5x5 SiC-SiC Interface GMC Table 5.4";
Vf = 0.28;% Fiber Volume Fraction
Vi  = 0.13;% Interface Volume Fraction
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 105;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);


Mat = 140;
Theory = 'GMC';
name = "26x26 GMC - Sun & Vaidya Carbon Epoxy";
Vf = 0.6;
Constits.Fiber = constitprops{13};
Constits.Matrix = constitprops{14};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 141;
Theory = 'HFGMC';
name = "26x26 HFGMC - Sun & Vaidya Carbon Epoxy";
Vf = 0.6;
Constits.Fiber = constitprops{13};
Constits.Matrix = constitprops{14};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 142;
name = "All theories glass-epoxy Vf = x";
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
% -- Call function (contained in this file) to get eff props across all Vf
RunOverVf(name, Constits, props, Mat);
props{Mat}.Quit = true; % -- Tells MicroAnalysis to quit after getting eff props

Mat = 143;
Theory = 'HFGMC';
name = "32x32 Glass-Epoxy Contour HFGMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 144;
Theory = 'HFGMC';
name = "100x100 Glass-Epoxy Contour HFGMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 99;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 145;
Theory = 'GMC';
name = "26x26 Glass-Epoxy Contour GMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 146;
Theory = 'GMC';
name = "68x120 Hex Glass-Epoxy Contour GMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 147;
Theory = 'HFGMC';
name = "68x120 Hex Glass-Epoxy Contour HFGMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 148;
Theory = 'GMC';
name = "79x79 Random (json) Glass-Epoxy Contour GMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 149;
Theory = 'HFGMC';
name = "79x79 Random (json) Glass-Epoxy Contour HFGMC";
Vf = 0.50;
Constits.Fiber = constitprops{8};
Constits.Matrix = constitprops{9};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 150;
Theory = 'HFGMC';
name = "136x236 Hex SiC-SiC HFGMC";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 151;
Theory = 'HFGMC';
name = "Grid convergence Hex Vf = 0.6 IM7-8552 HFGMC";
Vf = 0.6;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 152;
Theory = 'HFGMC';
name = "Random Vf = 0.5 IM7-8552 HFGMC";
Vf = 0.5;
Constits.Fiber = constitprops{1};
Constits.Matrix = constitprops{3};
RUCid = 300;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 153;  
Theory = 'GMC';
name = "26x26 Glass-Epoxy GMC";
Vf = 0.6;
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 154;  
Theory = 'HFGMC';
name = "26x26 Glass-Epoxy HFGMC";
Vf = 0.6;
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);


Mat = 171;
Theory = 'GMC';
name = "2x2 Glass-Epoxy GMC Vf = 0.55";
Vf = 0.55;
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
RUCid = 2;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 172;
Theory = 'HFGMC';
name = "notch in monolithic epoxy 25x25";
Vf = 0.2;
Constits.Fiber = constitprops{16};
Constits.Matrix = constitprops{15};
RUCid = 30;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 173;
Theory = 'HFGMC';
name = "notch in monolithic epoxy 51x51";
Vf = 0.2;
Constits.Fiber = constitprops{16};
Constits.Matrix = constitprops{15};
RUCid = 31;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 174;
Theory = 'HFGMC';
name = "notch in monolithic epoxy 101x101";
Vf = 0.2;
Constits.Fiber = constitprops{16};
Constits.Matrix = constitprops{15};
RUCid = 32;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 175;
Theory = 'HFGMC';
name = "notch in monolithic epoxy 201x201";
Vf = 0.2;
Constits.Fiber = constitprops{16};
Constits.Matrix = constitprops{15};
RUCid = 33;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 176;
Theory = 'MT';
name = "Section 7.3 - MT with MT props";
Vf = 0.55;
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 177;
Theory = 'HFGMC';
name = "Section 7.3 - HF with MT props";
Vf = 0.55;
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 178;
Theory = 'MT';
name = "Section 7.3 - MT with HF props";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 179;
Theory = 'HFGMC';
name = "Section 7.3 - HF with HF props";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 180;
Theory = 'GMC';
name = "Uni PMC - 32x32 GMC";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 181;
Theory = 'HFGMC';
name = "Uni PMC - 32x32 HFGMC";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 182;
Theory = 'HFGMC';
name = "Uni PMC - Random from json HFGMC";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 183;
Theory = 'HFGMC';
name = "Random HFGMC";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 300;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 184;
Theory = 'HFGMC';
name = "Uni PMC - Random from json HFGMC";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 185;
Theory = 'HFGMC';
name = "HF 7x7";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 7;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 186;
Theory = 'HFGMC';
name = "HF 32x32";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
RUCid = 26;
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);

Mat = 187;
Theory = 'MT';
name = "MT with HF props";
Vf = 0.55;
Constits.Fiber = constitprops{18};
Constits.Matrix = constitprops{15};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 188;
Theory = 'MT';
name = "MT with MT props";
Vf = 0.55;
Constits.Fiber = constitprops{2};
Constits.Matrix = constitprops{4};
[props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});

Mat = 189;
Theory = 'HFGMC';
name = "HFGMC 32x58 Hex CMC";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 190;
Theory = 'HFGMC';
name = "HFGMC CMC Random-A from json with int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 191;
Theory = 'HFGMC';
name = "HFGMC CMC Random-B from json with int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 192;
Theory = 'HFGMC';
name = "HFGMC CMC Random-C from json with int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 193;
Theory = 'HFGMC';
name = "HFGMC CMC Random-A from json, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 194;
Theory = 'HFGMC';
name = "HFGMC CMC Random-B from json, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 195;
Theory = 'HFGMC';
name = "HFGMC CMC Random-C from json, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 196;
Theory = 'HFGMC';
name = "HFGMC 32x58 Hex CMC, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 200;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 197;
Theory = 'HFGMC';
name = "HFGMC 13x13 square pack CMC with int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 113;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 198;
Theory = 'HFGMC';
name = "HFGMC CMC 4-Fiber from json with int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 199;
Theory = 'HFGMC';
name = "HFGMC CMC 8-Fiber from json with int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{12};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 200;
Theory = 'HFGMC';
name = "HFGMC 13x13 square pack CMC, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 113;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 201;
Theory = 'HFGMC';
name = "HFGMC CMC 4-Fiber from json, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);

Mat = 202;
Theory = 'HFGMC';
name = "HFGMC CMC 8-Fiber from json, NO int failure";
Vf = 0.28;
Vi = 0.13;
Constits.Fiber = constitprops{10};
Constits.Matrix = constitprops{11};
Constits.Interface = constitprops{17};
RUCid = 1000;
[props{Mat}] = RunMicro(Theory, name, struct('Vf',Vf,'Vi',Vi), Constits, props{Mat}, RUCid);


%========================================================================
% -- End of composite material definitions
%========================================================================

end

%========================================================================
% -- End composite material definitions
%========================================================================

%--------------------------------------------------------------
%--------------------------------------------------------------

function RunOverVf(name, Constits, props, Mat)

% -- Run Voigt, Reuss, MT, MOC, and MOCu micromechanics theories over all
%    Vf for the specifed Constits, write out eff prop results to file

Fact1 = 1000; % -- Factor - MPa to GPa
Fact2 = 1E6;  % -- Factor - /C to 1.E-6/C

VfResults = cell(5);

IncludeHFGMC = false;
%IncludeHFGMC = true;

if IncludeHFGMC 
   NumTheories = 9;
else
   NumTheories = 8;    
end

% -- Loop over theory
for k = 1:NumTheories
    switch k
        case 1
            Theory = 'Voigt';
            Name(k) = "Voigt       ";
        case 2
            Theory = 'Reuss';
            Name(k) = "Reuss       ";
        case 3
            Theory = 'MT';
            Name(k) = "MT          ";
        case 4
            Theory = 'MOC';
            Name(k) = "MOC         ";
        case 5
            Theory = 'MOCu';
            Name(k) = "MOCu        ";
        case 6
            Theory = 'GMC';
            RUCid = 2;
            Name(k) = "GMC-2x2     ";
        case 7
            Theory = 'GMC';
            RUCid = 7;
            Name(k) = "GMC-7x7     ";
        case 8
            Theory = 'GMC';
            RUCid = 26;
            Name(k) = "GMC-26x26   ";
        case 9
            Theory = 'HFGMC';
            RUCid = 26;
            Name(k) = "HFGMC-26x26 ";
    end

    N_Vf = 100;
    for i = 1: N_Vf + 1
        Vf = (i - 1)/N_Vf;
        if exist('RUCid','var')
            [props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat}, RUCid);
        else
            [props{Mat}] = RunMicro(Theory, name, Vf, Constits, props{Mat});
        end
        if ~props{Mat}.used
            return;
        end
        disp(['Vf = ', num2str(Vf)]);
        VfResults{k}.Vf(i) = Vf;
        VfResults{k}.EffProps(1,i) = props{Mat}.E1/Fact1;
        VfResults{k}.EffProps(2,i) = props{Mat}.E2/Fact1;
        VfResults{k}.EffProps(3,i) = props{Mat}.E3/Fact1;
        VfResults{k}.EffProps(4,i) = props{Mat}.G23/Fact1;
        VfResults{k}.EffProps(5,i) = props{Mat}.G13/Fact1;
        VfResults{k}.EffProps(6,i) = props{Mat}.G12/Fact1;
        VfResults{k}.EffProps(7,i) = props{Mat}.v12;
        VfResults{k}.EffProps(8,i) = props{Mat}.v13;
        VfResults{k}.EffProps(9,i) = props{Mat}.v23;
        VfResults{k}.EffProps(10,i) = props{Mat}.a1*Fact2;
        VfResults{k}.EffProps(11,i) = props{Mat}.a2*Fact2;
        VfResults{k}.EffProps(12,i) = props{Mat}.a3*Fact2;
    end
end

if ~IncludeHFGMC
    VfResults{9}.Vf = VfResults{8}.Vf;
    VfResults{9}.EffProps(1:12, 1:N_Vf + 1) = 0;
    Name(9) = "HFGMC Placeholder Ch6";
end

% -- Print to file for Excel plotting 
[~,~,~] = mkdir('Output');
tttt = datetime(datetime,'Format','yyyy-MMM-dd HH.mm.ss');
OutFile = ['Output/Props vs Vf - ',char(name),char(tttt),'.txt'];
fid = fopen(OutFile,'wt');
Fmet = '%E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \n';
Fmst = '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n';
Fmnt = '\n %s \n';
Props = ["E1", "E2", "E3", "G23", "G13", "G12", "V12", "V13", "V23", ...
         "alpha1", "alpha2", "alpha3"];

for j = 1:12
    TSProp = ['---- ', char(Props(j)),' ----'];
    fprintf(fid, Fmnt, TSProp);
    fprintf(fid, Fmst, ["Vf         ", Name]);
    for i = 1: N_Vf + 1
       fprintf(fid, Fmet, VfResults{1}.Vf(i), VfResults{1}.EffProps(j,i), ...
               VfResults{2}.EffProps(j,i), VfResults{3}.EffProps(j,i), ...
               VfResults{4}.EffProps(j,i), VfResults{5}.EffProps(j,i), ...
               VfResults{6}.EffProps(j,i), VfResults{7}.EffProps(j,i), ...
               VfResults{8}.EffProps(j,i), VfResults{9}.EffProps(j,i));
    end
end
   
fclose(fid);

end
