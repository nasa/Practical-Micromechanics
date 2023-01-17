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
% Purpose: This script runs stand-alone micromechanics progressive damage simulations  
%          treating the composite as a material point.  The composite material is 
%          defined in GetEffProps.m and the loading and problem setup is specified in 
%          MicroProblemDef.m
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Clear memory and close files
clear;
close all;
fclose('all');
clc;

% -- Add needed function locations to the path
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');
addpath('Functions/Damage');

%-----------------------------------------------------------------
% 1) Define Micromechanics Problems
%-----------------------------------------------------------------
[NProblems, OutInfo, MicroMat, Loads] = MicroProblemDef();

% -- Initializations
effprops = cell(1,300);
effprops_init = cell(1,300);
Results = cell(1,NProblems);

% -- Perform checks on problem definitions
for NP = 1: NProblems
    if MicroMat{NP} > 1
        
    else
       error(['No material defined for problem #', char(num2str(NP))]); 
    end
    
    if ~isfield(Loads{NP}, 'DT')
       Loads{NP}.DT = 0;
    elseif Loads{NP}.DT ~= 0
       error(['Thermal problem not implemented for progressive damage, DT ~= 0 in problem #', char(num2str(NP))]); 
    end

    if OutInfo.Format == "doc" 
        disp(' ---> Word output not supported for progressive damage, switching to txt');
        OutInfo.Format = "txt";
    end
    
    % -- Check for missing problem name
    if (ismissing(OutInfo.Name(NP)))
        OutInfo.Name(NP) = string(['Problem ', char(num2str(NP))]);
    end
    
end

%-----------------------------------------------------------------
% 2) Get constituent properties
%-----------------------------------------------------------------
[constitprops] = GetConstitProps;

%-----------------------------------------------------------------
% 3) Get inital effective properties from micromechanics
%-----------------------------------------------------------------

% -- Determine which mats are used in any problem
for NP = 1: NProblems
    mat = MicroMat{NP};
    effprops{mat}.used = true;
end

[effprops] = GetEffProps(constitprops, effprops);


% -- Loop through problems
for NP = 1: NProblems
    
    % -- Clear storage of stress/strain history for next problem
    if NP > 1
        clear strainhist;
        clear stresshist;
        clear fullhist;
    end
    
    OutInfo.PlottedInit = false; % -- Flag for plotting damage init location
    
    mat = MicroMat{NP};
    
    % -- Echo problem info to command window
    disp(['Micro Problem #',num2str(NP),' - ', char(OutInfo.Name(NP))]);
    disp(['   Material Number ',num2str(mat)]);    

    % -- Check that the problem's ply material has been defined
    if ~isfield(effprops{mat}, 'name')
        error(strcat('effective material #', num2str(mat), ' undefined ... check GetEffProps'));
    end
    
    % -- Get new RUC and new effprops if multiple problems and random RUC
    if NP > 1
        if isfield(effprops{mat}, 'RUC')
            RUCid = effprops{mat}.RUC.id;
            if RUCid == 300
                Theory = char(effprops{mat}.micro);
                if isfield(effprops{mat}, 'Vi')
                    Constits.Interface = effprops{mat}.Interface;
                else
                    effprops{mat}.Vi = 0;
                end
                Constits.Fiber = effprops{mat}.Fiber;
                Constits.Matrix = effprops{mat}.Matrix;
                Vol = struct('Vf', effprops{mat}.Vf, 'Vi', effprops{mat}.Vi);
                [effprops{mat}] = RunMicro(Theory, effprops{mat}.name, Vol, Constits, effprops{mat}, RUCid);
            end
        end
    end

    % -- Store initial eff props and RUC -- this will get changed when damage occurs
    effprops_init{mat} = effprops{mat};
    
    % -- Specifed loads are final values for incremental loading
    Loads{NP}.Final = Loads{NP}.Value;
    Loads{NP}.DTFinal = Loads{NP}.DT;
    
% -- Prog dam not enabled for MOC or MOCu -- these are implemented for fiber/matrix not 
%    different materials per matrix subcell
    if effprops{mat}.micro == "MOC" || effprops{mat}.micro == "MOCu"
        error('        Progressive damage not implemented for MOC and MOCu')
    end
    
%-----------------------------------------------------------------
% 4) Create Damaged Constituents
%-----------------------------------------------------------------
    DamageFactor = 0.0001; % -- Stiffness reduction factor for damaged material
    [effprops{mat}.FiberDam] = DamageConstit(effprops{mat}.Fiber, DamageFactor);
    [effprops{mat}.MatrixDam] = DamageConstit(effprops{mat}.Matrix, DamageFactor);
    if isfield(effprops{mat},'Interface')
        [effprops{mat}.InterfaceDam] = DamageConstit(effprops{mat}.Interface, DamageFactor);
    end
    
    % -- Default to 100 increments if NINC not specified
    if ~isfield(Loads{NP}, 'NINC')
        Loads{NP}.NINC = 100;
    end

    % -------------------------
    % -- Loading increment loop
    % -------------------------
    INC = -1;
    INCcount = 0;
    while INC < Loads{NP}.NINC
    
        INC = INC + 1;
        INCcount = INCcount + 1;
        OutInfo.INC = INC;
        
        % -- Determine current load level
        pct = INC/Loads{NP}.NINC;
        Loads{NP}.Value = pct*Loads{NP}.Final;
        Loads{NP}.DT = pct*Loads{NP}.DTFinal;   
        
        % -- Display increment progress 
        disp(' ');
        disp('*****************************************************');
        disp(['* INC = ',num2str(INC), '  out of ', num2str(Loads{NP}.NINC)]);
        disp('*****************************************************');
        
        % -----------------
        % -- Iteration Loop
        % -----------------
        NITER = 1000000000; % -- Max number of iterations
        for ITER = 1:NITER
            
            % -- Display iteration progress 
            disp(' ');
            disp(['----- ITER = ', char(num2str(ITER)), ' -----']);
            OutInfo.ITER = ITER;
            
            %-----------------------------------------------------------------
            % 5) Solve loading and calculate local fields for micromechanics
            %-----------------------------------------------------------------
            epsth = Loads{NP}.DT * [effprops{mat}.a1; effprops{mat}.a2; effprops{mat}.a3; 0; 0; 0];
            B = -effprops{mat}.Cstar*epsth;
            [SG, FullGlobalStrain] = SolveLoading(6, effprops{mat}.Cstar, B, Loads{NP});

            if ~exist('Results', 'var') 
                [Results{NP}] = MicroFields(FullGlobalStrain, Loads{NP}.DT, effprops{mat});
            elseif length(Results) < NP
                [Results{NP}] = MicroFields(FullGlobalStrain, Loads{NP}.DT, effprops{mat});           
            else
                [Results{NP}] = MicroFields(FullGlobalStrain, Loads{NP}.DT, effprops{mat}, Results{NP});
            end

            Stress = effprops{mat}.Cstar*(FullGlobalStrain - epsth);
            
            %------------------------------------------------------------------
            % 6) Write initial micromechanics properties and write/plot results
            %------------------------------------------------------------------
            
            if INC == 0
                effprops{mat}.Mat = mat;
                [OutInfo] = OutputMicro(OutInfo, NP, effprops{mat}, Loads{NP}, SG, FullGlobalStrain, epsth, true);

                % -- Inverse of Cstar
                CstarInv = inv(effprops{mat}.Cstar);

                % -- Store initial trace(inv Cstar) - used as global compliance measure
                InitTraceCI = CstarInv(1,1) + CstarInv(2,2) + CstarInv(3,3);
                TraceCI = InitTraceCI;
            end

            % -- Loads{NP}.PlotComponent can be specified in MicroProblemDef.m, or 
            %    will plot applied direction if only 1 component applied
            if isfield(Loads{NP}, 'PlotComponent')
                PlotComponent = Loads{NP}.PlotComponent;
            else
                PlotComponent = 2;
                for i = 1: 6
                    if Loads{NP}.Final(i) ~= 0
                        PlotComponent = i;
                        break;
                    end
                end
            end
            
            % -- Store global stress and strain history
            strainhist(INCcount) = FullGlobalStrain(PlotComponent);
            stresshist(INCcount) = Stress(PlotComponent);
            fullhist(INCcount).Stress = Stress;
            fullhist(INCcount).Strain = FullGlobalStrain;
            if INC == 0
                damfig = 0;
                Results{NP}.MoS.NewFailure = false;
            end
            
            % -- Function (included below) to write damage results
            [damfig, strainhist, stresshist, OutInfo] = DamageWrite(NP, OutInfo, effprops{mat}, Results{NP}, Stress, FullGlobalStrain, INC, ...
                                                       ITER, strainhist, stresshist, fullhist, damfig);            
            
            
            % -- Stop execution for if TraceCI exceeds termination factor
            if isfield(Loads{NP}, 'TerminationFactor')
                factor = Loads{NP}.TerminationFactor;
            else
                factor = 50000; % -- Defaults to very high factor
            end
            if TraceCI/InitTraceCI > factor
                INC = Loads{NP}.NINC + 1;
                break
            end

            %-----------------------------------------------------------------
            % 7) Calculate Margins
            %-----------------------------------------------------------------
            % -- Check for turning off some failure criteria (Tsai-Wu only by default, 0 = off)
            if (~isfield(Loads{NP},'CriteriaOn'))
                Loads{NP}.CriteriaOn = [0,0,0,1];
            end
            [Results{NP}] = CalcMicroMargins(effprops{mat}, Loads{NP}, Results{NP}, true); 

            % -- Initialize failure flags
            if INC == 0
                if Results{NP}.Type == "FM"
                    Results{NP}.MoS.MicroF.Failed = false;
                    Results{NP}.MoS.MicroM.Failed = false;
                else
                    for b = 1: effprops{mat}.RUC.NB
                        for g = 1: effprops{mat}.RUC.NG
                            Results{NP}.MoS.Micro(b,g).Failed = false;
                        end 
                    end
                end
            end

            %-----------------------------------------------------------------------------
            % 8) Reduce stiffness of any damaged subcell (those within a tolerance of the  
            %    min negative margin)
            %-----------------------------------------------------------------------------
            [effprops{mat}, Results{NP}.MoS] = DamageMicro(effprops{mat}, Results{NP}.MoS, OutInfo);
           
            % -- Find min MoS to skip (jump) increments
            CheckForSkip = true;
            if ITER ~= 1
                CheckForSkip = false;
            elseif isfield(Results{NP}.MoS, 'NewFailure')
                if Results{NP}.MoS.NewFailure
                    CheckForSkip = false;
                end
            end
            
            if CheckForSkip
                MinMoS = 99999;
                if Results{NP}.Type == "FM"
                    if isfield(Results{NP}.MoS.MicroF, 'MinMoS')
                        MinMoS = min(MinMoS, Results{NP}.MoS.MicroF.MinMoS);
                    end
                    if isfield(Results{NP}.MoS.MicroM, 'MinMoS')
                        MinMoS = min(MinMoS, Results{NP}.MoS.MicroM.MinMoS);
                    end
                else
                    for b = 1: effprops{mat}.RUC.NB
                        for g = 1: effprops{mat}.RUC.NG

                            if ~isfield(Results{NP}.MoS.Micro(b,g), 'MinMoS')
                                continue;
                            end
                            
                            if isfield(Results{NP}.MoS.Micro(b,g), 'Failed') 
                                if Results{NP}.MoS.Micro(b,g).Failed
                                    continue;
                                else
                                    MinMoS = min(MinMoS, Results{NP}.MoS.Micro(b,g).MinMoS);
                                end
                            else
                                MinMoS = min(MinMoS, Results{NP}.MoS.Micro(b,g).MinMoS);
                            end
                        end 
                    end
                end

                disp(['  MinMoS = ',char(num2str(MinMoS))]);
                
                % -- Step out to just before increment where next damage should occur
                %    if MinMoS is positive (no additional failures)
                if MinMoS > 0 && MinMoS < 99999 && INC < Loads{NP}.NINC
                    INC_JUMP = round( INC*( (MinMoS + 1) - 1) );
                    disp(['  INC_JUMP = ', char(num2str(INC_JUMP))]);
                    if INC_JUMP > 2
                        INC = INC + INC_JUMP - 2;
                        % -- Jump to final increment if will no more damage will occur
                        if INC > Loads{NP}.NINC - 1
                            INC = Loads{NP}.NINC - 1;
                        end
                        break;
                    end
                elseif MinMoS == 99999 && INC > 10 && INC < Loads{NP}.NINC
                    INC = Loads{NP}.NINC - 1;
                    break;
                end
             end

            % -- Inverse of Cstar
            CstarInv = inv(effprops{mat}.Cstar);

            % -- Check for large stiffness drop (Cinv increase)
            TraceCI = CstarInv(1,1) + CstarInv(2,2) + CstarInv(3,3);

            disp(['  --> Cinv Trace Factor = ',char(num2str(TraceCI/InitTraceCI))]);
        
            % -- Exit iteration loop if no new failure
            if ~Results{NP}.MoS.NewFailure
                break
            end
            
        end % -- End Iteration Loop
        
    end % -- End Increment Loop
        
    %-------------------------------------------------------------------------
    % 9) Plot final micro fields and write RUC  (only for RUCid = 300 or 1000)
    %-------------------------------------------------------------------------
    if (~isfield(OutInfo,'Format'))
        OutInfo.Format = "txt";
    end

    if OutInfo.MakePlots 
        PlotMicroFields(OutInfo, effprops{mat}, Results{NP});
    end

    WriteRUC(OutInfo, OutInfo.Name(NP), effprops{mat})
    
    
    % -- Reset effprops and RUC to init values 
    %   (in case this material is used in another problem)
    effprops{mat} = effprops_init{mat};
    
    
    disp(['  *** Problem ',char(num2str(NP)),' Completed ***'])
    disp(' ');

    figh = findall(0,'type','figure');
    other_figures = setdiff(figh, damfig);
    delete(other_figures)
    
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [damfig, strainhist, stresshist, OutInfo] = ...
    DamageWrite(NP, OutInfo, props, Results, Stress, FullGlobalStrain, INC, ...
                ITER, strainhist, stresshist, fullhist, damfig)

% -- Write and plot progressive damage results as damage progresses

MoS = Results.MoS;
if isfield(Results, 'MicroFields')
    MicroFields = Results.MicroFields;
end

% -- Write stress-strain for every iteration (including non-converged states)
OutInfo.DamDir = ['Output/DamSnapShots - ',char(OutInfo.Name(NP)),' ',OutInfo.datetime,'/'];
[~,~,~] = mkdir(OutInfo.DamDir);
filename = [OutInfo.DamDir, 'ITER stress vs strain - ', char(OutInfo.Name(NP)),'.txt']; 
fid = fopen(filename, 'a');

if INC == 0
    Fmt1 = '\r %s';
    text = 'INC ITER strain11 stress11 strain22 stress22 strain33 stress33 strain23 stress23 strain13 stress13 strain12 stress12';
    style='Normal';
    WriteText('txt', fid, Fmt1, 0, text, style, [0,0])
end

Fmt1 = '\r %s';
text = [num2str(INC), ' ', num2str(ITER), ' ', ...
        num2str(FullGlobalStrain(1)),' ',num2str(Stress(1)),' ',num2str(FullGlobalStrain(2)),' ',num2str(Stress(2)),' ', ...
        num2str(FullGlobalStrain(3)),' ',num2str(Stress(3)),' ',num2str(FullGlobalStrain(4)),' ',num2str(Stress(4)),' ',...
        num2str(FullGlobalStrain(5)),' ',num2str(Stress(5)),' ',num2str(FullGlobalStrain(6)),' ',num2str(Stress(6))];
style='Normal';
WriteText('txt', fid, Fmt1, 0, text, style, [0,0])

fclose(fid);

% -- Write stress-strain for only last iteration (no non-converged states)
OutInfo.DamDir = ['Output/DamSnapShots - ',char(OutInfo.Name(NP)),' ',OutInfo.datetime,'/'];
[~,~,~] = mkdir(OutInfo.DamDir);
filename = [OutInfo.DamDir, 'stress vs strain - ', char(OutInfo.Name(NP)),'.txt']; 
fid = fopen(filename,'w');

Fmt1 = '\r %s';
text = 'strain11 stress11 strain22 stress22 strain33 stress33 strain23 stress23 strain13 stress13 strain12 stress12';
style='Normal';
WriteText('txt', fid, Fmt1, 0, text, style, [0,0])

Fmt1 = '\r %s';
style='Normal';
for I = 1: length(fullhist)
    text = [num2str(fullhist(I).Strain(1)),' ',num2str(fullhist(I).Stress(1)),' ',num2str(fullhist(I).Strain(2)),' ',num2str(fullhist(I).Stress(2)),' ', ...
            num2str(fullhist(I).Strain(3)),' ',num2str(fullhist(I).Stress(3)),' ',num2str(fullhist(I).Strain(4)),' ',num2str(fullhist(I).Stress(4)),' ',...
            num2str(fullhist(I).Strain(5)),' ',num2str(fullhist(I).Stress(5)),' ',num2str(fullhist(I).Strain(6)),' ',num2str(fullhist(I).Stress(6))];
    WriteText('txt', fid, Fmt1, 0, text, style, [0,0])
end

fclose(fid);


% -- Write eff props per increment and ITER (including non-converged states)
filename = [OutInfo.DamDir, 'ITER props vs INC - ', char(OutInfo.Name(NP)),'.txt']; 
fid = fopen(filename, 'a');
Fmt1 = '\r %s';
if INC == 0
    Fmt1 = '\r %s';
    text = 'INC ITER E11 E22 E33 v12 v13 v23 G12 G13 G23 alpha1 alpha2 alpha3';
    style='Normal';
    WriteText('txt', fid, Fmt1, 0, text, style, [0,0])
else
    text = [num2str(INC),' ',num2str(ITER),' ',num2str(props.E1),' ',num2str(props.E2),' ',num2str(props.E3),' ',num2str(props.v12),' ', ...
            num2str(props.v13),' ',num2str(props.v23),' ',num2str(props.G12),' ',num2str(props.G13),' ',...
            num2str(props.G23),' ',num2str(props.a1),' ',num2str(props.a2),' ',num2str(props.a3)];
    style='Normal';
    WriteText('txt', fid, Fmt1, 0, text, style, [0,0])
end
fclose(fid);

% -- Plot stress vs. strain
if INC == 0
    damfig = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
else
    figure(damfig);
    subplot(1,3,1)
    plot(strainhist,stresshist,'-o')
    title('stress vs. strain');
    xlabel('strain'); 
    ylabel('stress (MPa)','rotation',90);
    drawnow
end

% -- Display and save stress-strain plot, local von Mises stress and damage

% -- Only update/write plot when new failure happens
if (MoS.NewFailure || INC == 1) && (props.micro == "GMC" || props.micro == "HFGMC")
% -- Plot/write every increment
%if INC > 0 && (props.micro == "GMC" || props.micro == "HFGMC")

    mats = zeros(2*props.RUC.NB, 2*props.RUC.NG);
    x2 = zeros(1, 2*props.RUC.NB);
    x3 = zeros(1, 2*props.RUC.NG);
    Svm = zeros(2*props.RUC.NB, 2*props.RUC.NG);
    
    % -- Plot vm stress and failed subcells
    for j = 1:2*props.RUC.NG
        g = round(j/2);
        for i = 1:2*props.RUC.NB
            b = round(i/2);
            switch props.RUC.matsCh(b,g)
                case 'F'
                    mats(i,j) = 1;
                case 'M'
                    mats(i,j) = 2;
                case 'I'
                    mats(i,j) = 1.6;
                otherwise
                    mats(i,j) = 3;
            end
        end 
    end 

    % place an x2-point at bottom and top of each subcell
    h = 0;
    i = 1;
    for b = 1:props.RUC.NB
        x2(i) = h;
        h = h + props.RUC.h(b);
        i = i + 1;
        x2(i) = h;
        i = i + 1;
    end

    % place an x3-point at left and right of each subcell
    l = 0;
    i = 1;
    for g = 1:props.RUC.NG
        x3(i) = l;
        l = l + props.RUC.l(g);
        i = i + 1;
        x3(i) = l;
        i = i + 1;
    end
    
   % -- Create grid
    [X,Y] = meshgrid(x3,x2);

   % -- Calc von Mises stress
    for j = 1:2*props.RUC.NG
        g = round(j/2);
        for i = 1:2*props.RUC.NB
            b = round(i/2);
            stress = MicroFields(b,g).stress;
            Svm(i,j) = sqrt( (stress(1) - stress(2))^2 + ...
                             (stress(2) - stress(3))^2 + ...
                             (stress(1) - stress(3))^2 + ...
                              6*stress(4)^2 + ...
                              6*stress(5)^2 + ...
                              6*stress(6)^2 ) / sqrt(2);
        end 
    end 

    % --- Plot von Mises ---
    subplot(1,3,2)
    colormap(jet);
    pcolor( X, Y, Svm), shading interp;
    %pcolor( X, Y, Svm), shading faceted;
    fontsize = 12;
    c=colorbar;
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title(['von Mises stress, INC = ',char(num2str(INC)),' ITER = ',char(num2str(ITER))],'FontSize',fontsize);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off; 
    drawnow

    % --- Plot damage ---
    subplot(1,3,3)
    colormap(jet);
    %pcolor( X, Y, mats), shading interp;
    pcolor( X, Y, mats), shading faceted;
    caxis([1 3]);
    fontsize = 12;
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title(['Subcell Failure, INC = ',char(num2str(INC)),' ITER = ',char(num2str(ITER))],'FontSize',fontsize);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off; 
    drawnow

    % -- Save damage snapshot figure
    if OutInfo.Format == "doc"
         WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
    elseif OutInfo.Format == "txt"
        if ~OutInfo.PlottedInit
            if ITER == 2
                Plotname = [OutInfo.DamDir,'Initiation Location Damage Snapshot  - ', 'INC = ', char(num2str(INC)), ' ITER = ', char(num2str(ITER)),'.bmp'];
                saveas(gcf,Plotname);
                OutInfo.PlottedInit = true;
            end
        end
        % -- This will overwrite the plot each ITER for a given INC, so left with final
        Plotname = [OutInfo.DamDir,'Damage Snapshot  - ', 'INC = ', char(num2str(INC)),' FINAL','.bmp'];
        saveas(gcf,Plotname);
        
        % -- This will write a plot for every ITER -- so can produce a lot of data
%         Plotname = [OutInfo.DamDir,'Damage Snapshot  - ', 'INC = ', char(num2str(INC)), 'ITER = ', char(num2str(ITER)),'.bmp'];
%         saveas(gcf,Plotname);
    end      

end

end
