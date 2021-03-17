function [fvr,ivr,RUC] = RUCRand(Fvr,Radf,Nfibers,Tint,touching,MaxTries,Plot)
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

% Program to generate random placement of circular fiber RUC
% in a nx x ny Region...
% Developed by Pappu L.N. Murthy, NASA GRC
% Ref: Murthy, P.L.N., Bednarcyk, B.A., and Mital, S.K. (2017) “A Compilation of MATLAB 
%      Scripts and Functions for MAC/GMC Analysis” NASA/TM-2017-219500.

Radi = Radf + Tint;
if (touching)
    dpad = 0;
else
    dpad = 1;
end

M=2;F=1;I=3;U=5;
e=0;
[BaseRUC]=BaseCmcRuc3(Radf,Radi,e);
[Dia,~]=size(BaseRUC);
Rad=Dia/2;
FA = length ( BaseRUC(BaseRUC==F));
IA = length ( BaseRUC(BaseRUC==I));

% Find Big RUC size
Side= round( sqrt( FA*Nfibers/Fvr));
MA = Side^2-Nfibers*(FA+IA);
fvr = FA*Nfibers/Side/Side;
ivr = IA*Nfibers/Side/Side;
MVR = MA/Side/Side;
nx=Side;
ny=Side;
%--------------------------------------------------------------------------
% dpad=0; % Fibers can touch
% dpad = 1, minimum of one subcell between the fibers.
IndForI=find(BaseRUC==F | BaseRUC==I);
% Mark all the nodes that will be unavailable after placing a fiber
ElimRUC=BaseCmcRuc3(2*Radf,2*(Radi+dpad),0.1);
ElimRUC(ElimRUC==F | ElimRUC==I)=U;
IndForU = find(ElimRUC==U);
%--------------------------------------------------------------------------
% All locations in ElimRUC are unavailable once a BaseRUC is placed in
% place at any location. The positions of the subcells change depending on
% the fiber center.
%--------------------------------------------------------------------------
[NsizeE,~]=size(ElimRUC);
Maxgen=MaxTries;
RUC=M*ones(nx,ny);
RUCU=RUC;

% Initialize variables
GenV=zeros(1,Maxgen);
NfibV=zeros(1,Maxgen);
Nf =zeros(Nfibers,2); % Fiber Centers
nfib = 0;

% Define the domain for random center
xmin = 1; xmax = nx ;
ymin = 1; ymax = ny ;

%--------------------------------------------------------------------------
% Make a list of all eligible points in the domain.
NsizeB=(xmax-xmin+1);
EligRefL =cell(NsizeB^2,1);
EligIndx = (1:NsizeB^2)';
X=xmin:xmax;
Y=ymin:ymax;
l=0;
for j = 1:NsizeB
    for i=1:NsizeB
        l=l+1;
        EligRefL {l}= [X(i),Y(j), RUCU(i,j)];
    end
end
EligibleL=EligRefL;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Igen=1;
while ( Igen <= MaxTries && nfib < Nfibers) % First While Loop checking Number of fibers placed

    % Generate random center;
    GenV(Igen)=Igen;
    % generate uniform random number
    % pick a random location;
    rng(1000*second(datetime));
    Rloc = randi(length(EligIndx),1);
    LocR = EligibleL{Rloc};
    
    Cx= LocR(1,1);Cy=LocR(1,2);
    FibLoc = [Cx , Cy ];
    
    Nf(nfib+1,:) = FibLoc(1,:); % Store Fiber Centers
    % show square fiber for now
    px=FibLoc(1,1)-Rad;py=FibLoc(1,2)-Rad;
    Irow = px:px+Dia-1;
    Jcol = py:py+Dia-1;
    
    % Check if TempRUC dimension will go outside of box EJP
    [Irow,Jcol]=CheckBDcross(ny,nx,Irow,Jcol);
    tmp3 = RUC(Irow,Jcol);
    tmp3(IndForI)=BaseRUC(IndForI);
    RUC(Irow,Jcol)=tmp3;
    % figure(1), plotRUC(RUC )
    
    % Add code for Unavailable locations here
    pxU=FibLoc(1,1)-NsizeE/2;pyU=FibLoc(1,2)-NsizeE/2;
    IrowU = pxU:pxU+NsizeE-1;
    JcolU = pyU:pyU+NsizeE-1;
    
    % Check if ElimRUC dimension will go outside of box
    [IrowU,JcolU]=CheckBDcross(ny,nx,IrowU,JcolU);
    ind = find(RUCU(IrowU,JcolU)==U);
    tmp = RUCU(IrowU,JcolU);
    tmp2 = ElimRUC;
    tmp2(ind) = tmp(ind);
    RUCU(IrowU,JcolU) = tmp2;
    
    % Increment the trial counter
    nfib=nfib+1;
    NfibV(Igen)=nfib;
    
    % Update Eligible Region by removing the fiber that is already
    % placed with in the domain
    [EligIndx,EligibleL] = PurgeAndUpdate(IrowU,NsizeE,JcolU,...
    NsizeB,EligIndx,EligRefL,ElimRUC);
    TF=isempty(EligIndx);
    
    if (TF && nfib < Nfibers)
    % Start the counter again
        Igen = Igen+1;
        nfib=0;
        EligIndx = (1:NsizeB^2)';
        EligibleL=EligRefL;
        RUC=M*ones(nx,ny);
        RUCU=RUC;
        Nf = zeros(Nfibers,2); % Fiber Centers
    end
    
end

% -------------------------------------------------------------------------
% Check Whether max generations reached...
if (nfib < Nfibers && Igen>MaxTries)
    Msg={'Max. Gens Reached','Failed to Place all Fibers',...
    'Decrease Fvr','Or increase Radf and try again'};
    warndlg(Msg,'Random RUC Gen Failed');
end
% -------------------------------------------------------------------------

Generation=GenV(1:Igen)';
N_of_fibers=NfibV(1:Igen)';
Ta= table(Generation,N_of_fibers);
disp('Random RUC Iterations');
disp(Ta);
%Lbase=nx;Hbase=ny;
%[BigRUC,BigL,BigH] = RepeatRUC (Lbase,Hbase,RUC,nrep);
FibCenters = [ Nf(:,1)-1 Nf(:,2)];

H = ones(1,nx);
L = ones(1,ny);

% -- plot
if (Plot)
    figure
    [m,n]=size(RUC); Nmat = length(unique(RUC));
    px=0;py=0;
    switch Nmat
        case 1
            cmap=[0 0.5 1];
        case 2
            cmap=[0 0.5 1;1 1 0.5];
        case 3
            cmap= [0 0.5 1;1 1 0.5;0 0 0];
        case 4
            cmap= [0 0.5 1;1 1 0.5;0 0 0;0 1 1];
    end
    for i = 1:m
        for j=1:n
            rectangle('Position',[px py L(j) H(i)],'FaceColor',...
            cmap(RUC( i ,j),:));
            px=px+L(j);
        end
        px=0;
        py=py+H(i);
    end
    xlim([0 sum(L )]);ylim([0 sum(H )]);
    daspect([1,1,1]);
    txt = [num2str(m), 'x' , num2str(n),' Random RUC, Vf = ', num2str(fvr), ', Vi = ', num2str(ivr)];
    title(txt);
    axis off;
    drawnow
%     dbox = warndlg('Click OK when finished inspecting RUC','RUC Inspection');
%     uiwait(dbox);
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function[EligIndx,EligibleL] = PurgeAndUpdate(Irow,NsizeS,Jcol,...
NsizeB,EligIndx,EligRef,ElimRUC)
% Update Eligible Region by removing the locations occupied by the fibers
% that are already placed with in the domain
l=0;
for j = 1:NsizeS
    for i= 1:NsizeS
        if(ElimRUC(i,j)==5)
            l=l+1;
            SmallL(l)=Irow(i) + (Jcol(j)-1)*NsizeB ;
        end
    end
end

% Purging operation
EligIndx=setdiff(EligIndx,SmallL);
EligibleL=EligRef(EligIndx);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function[base]=BaseCmcRuc3(Radf,Radi,e)
% Program that generates a Square packing CMC Ruc with interface.
% Fiber radius and Interface thickness are user definable.
% Developed by Pappu L.N. Murthy, Date June, 30, 2015;
%
M=2;F=1;I=3;
nx=round(Radf); ny=round(Radi);
if(ny>nx)
    nx=ny;
end
base = M*ones(nx ,ny ); % Square RUC

% define center of RUC
cx=nx;cy=0;
[X,Y]= meshgrid((1:nx)-0.5 , (1:ny)-0.5 ) ;
Dist = sqrt( (cx-X).^2 + (cy-Y).^2 )-e;
base(Dist<=Radf)=F;
base(Dist>Radf & Dist<=Radi) = I;
base= [flipud([ base,fliplr(base)]);[ base,fliplr(base)]];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function[Irow,Jcol]=CheckBDcross(ny,nx,Irow,Jcol)
% Check if TempRUC dimension will go outside of box EJP
% This function checks each fiber block to see if it crosses the boundary
% and if it crosses the boundary it tries to place the remaining fiber on
% the opposite side to preserve periodicity of the RUC

for i=1:length(Irow)
    if Irow(i)>ny
        Irow(i)=Irow(i)-ny;
    elseif Irow(i)<1
        Irow(i)=Irow(i)+ny;
    end
end

for j=1:length(Jcol)
    if Jcol(j)>nx
        Jcol(j)=Jcol(j)-nx;
    elseif Jcol(j)<1
        Jcol(j)=Jcol(j)+nx;
    end
end