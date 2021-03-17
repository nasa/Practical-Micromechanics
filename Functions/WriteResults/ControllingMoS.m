function [MoS] = ControllingMoS(Type, props, MoS)
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
% Purpose: Determine controlling failure criterion and store related info
% Input:
% - Type: Average constituent "FM", or subcell "bg" based fields
% - props: Struct containing composite material and constituent properties 
% - MoS: Struct containing MoS results and info 
% Output: 
% - MoS: Updated struct containing MoS results and info
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MoS.MinMoS = 99999;

if Type == "FM" % -- Average constituent fields
    MoS.MicroF.MinMos = 99999;
    MoS.MicroM.MinMos = 99999;
    for j = 1:4
        if MoS.MicroF.Crit{j}.MinMoS < MoS.MicroF.MinMos
           MoS.MicroF.MinMos = MoS.MicroF.Crit{j}.MinMoS;
           MoS.MicroF.ConrollingNum = j;
           MoS.MicroF.Conrolling = MoS.MicroF.Crit{j}.name;
           MoS.MicroF.ConrollingComp = MoS.MicroF.Crit{j}.Controlling;
        end
        if MoS.MicroM.Crit{j}.MinMoS < MoS.MicroM.MinMos
           MoS.MicroM.MinMos = MoS.MicroM.Crit{j}.MinMoS;
           MoS.MicroM.ConrollingNum = j;
           MoS.MicroM.Conrolling = MoS.MicroM.Crit{j}.name;
           MoS.MicroM.ConrollingComp = MoS.MicroM.Crit{j}.Controlling;
        end
    end
    if MoS.MicroF.MinMos < MoS.MinMoS
        MoS.MinMoS = MoS.MicroF.MinMos;
        MoS.Controlling = 'Fiber';
        MoS.ControllingType = 'FM';
        MoS.ControllingComp = MoS.MicroF.ConrollingComp;
    end
    if MoS.MicroM.MinMos < MoS.MinMoS
        MoS.MinMoS = MoS.MicroM.MinMos;
        MoS.Controlling = 'Matrix';
        MoS.ControllingType = 'FM';
    end

else % -- Subcell fields
    for b = 1:props.RUC.NB
       for g = 1:props.RUC.NG
            MoS.Micro(b,g).MinMos = 99999;
            for j = 1:4
               if MoS.Micro(b,g).Crit{j}.MinMoS < MoS.Micro(b,g).MinMos
                   MoS.Micro(b,g).MinMos = MoS.Micro(b,g).Crit{j}.MinMoS;
                   MoS.Micro(b,g).ConrollingNum = j;
                   MoS.Micro(b,g).Conrolling = MoS.Micro(b,g).Crit{j}.name;
                   MoS.Micro(b,g).ConrollingComp = MoS.Micro(b,g).Crit{j}.Controlling;
               end
            end
            if MoS.Micro(b,g).MinMos < MoS.MinMoS
                MoS.MinMoS = MoS.Micro(b,g).MinMos;
                MoS.Controlling = [b, g];
                MoS.ControllingType = 'bg';
            end                         
       end
    end

end

end