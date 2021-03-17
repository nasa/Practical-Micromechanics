function [OutInfo] = OutputLam(OutInfo, Problem, Geometry, Loads, plyprops, LamResults)
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
% Purpose: Write output for laminate problems
% Input: 
% - OutInfo: Struct containing output information
% - Problem: Current problem number
% - Geometry: Struct containing laminate definition variables
% - Loads: Struct containing problem loading information
% - plyprops: Cell array containing ply material properties
% - LamResults: Struct containing laminate analysis results
% Output:
% - OutInfo: Updated struct containing output information
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Extract variables for convenience
LayerZ = Geometry.LayerZ;
Orient = Geometry.Orient;
N = Geometry.N;
DT = Loads.DT;
NM = LamResults.NM;
stress = LamResults.stress;
sigmaX = stress(1,:);
sigmaY = stress(2,:);
TauXY = stress(3,:);

% -- Set up output filenames and paths
[OutInfo] = SetupOutput(OutInfo, Problem);

% -- Write output to word or text
WriteLamOutput(OutInfo, Problem, Geometry, Loads, LamResults, plyprops);

% - Plot ply level stress through the thickness
if OutInfo.MakePlots && ~isfield(Loads, 'ang') % -- Turn off plotting for envelopes

    FS = 12 ;
    LW = 2;
    FW = 'bold';

    hold on
    figure(Problem)
    %plottools on;
    plot(sigmaX, LayerZ,'g-', 'LineWidth', LW) 
    %axis ij
    set(gca,'FontWeight',FW,'FontSize', FS);
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel('stress');
    ylabel('z');
    h = gca;
    h.YDir = 'reverse';
    Fmt = repmat('%.1f, ',1,N-1 );

    title(['Lamina Stresses: ', OutInfo.Name(Problem), ...
        sprintf('NM = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]',NM)...
        sprintf('delta T = [%.2f]',DT)...
        sprintf(['Ply Angles = [', Fmt,'%.1f]'],Orient),...
        sprintf(' ')],...
        'FontSize',2*FS/3);

    hold on 
    figure(Problem)
    plot(sigmaY,LayerZ,'k', 'LineWidth', LW) 
    %axis ij
    xlabel('stress');
    ylabel('z');

    figure(Problem)
    plot(TauXY,LayerZ,'r--', 'LineWidth', LW);
    %axis ij
    xlabel('stress');
    ylabel('z');
    box on;
    lgd=legend ('SigmaX','Sigma Y','TauXY');
    lgd.Location='eastoutside';
    lgd.NumColumns=1;
    lgd.FontSize=2*FS/3;
    legend boxoff;
    pause(2);

    % -- If using Word format add plot files to Word doc else store as jpg
    if OutInfo.Format == "doc"
        WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
    elseif OutInfo.Format == "txt"
        Plotname = [OutInfo.PlotFile,'.jpg'];
        saveas(gcf,Plotname);
    end

end
