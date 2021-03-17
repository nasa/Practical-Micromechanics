function [CritName, Comp, b, g] = WriteControlling(Format, fid, Fmt1, Fmt2, ...
                                                   word, MoS, props)
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
% Purpose: Write the overall controlling MoS info for the problem
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

style='Normal';

if string(MoS.ControllingType) == "FM"
    text=['Controlling Constituent is ', MoS.Controlling];
    if string(MoS.Controlling) == "Fiber"
        CritNum = MoS.MicroF.ConrollingNum;
        CritName = MoS.MicroF.Conrolling;
        Comp = MoS.MicroF.ConrollingComp;
    else
        CritNum = MoS.MicroM.ConrollingNum;
        CritName = MoS.MicroM.Conrolling;
        Comp = MoS.MicroM.ConrollingComp;
    end
    b = 0;
    g = 0;

elseif string(MoS.ControllingType) == "bg"
    b = MoS.Controlling(1);
    g = MoS.Controlling(2);
    text=['Controlling Subcell is (',num2str(b),',',num2str(g),')',' - ',props.RUC.matsCh(b,g)];
    CritNum = MoS.Micro(b,g).ConrollingNum;
    CritName = MoS.Micro(b,g).Controlling;
    Comp = MoS.Micro(b,g).ConrollingComp;
end
WriteText(Format, fid, Fmt1, word, text, style, [1,0])

if CritNum < 3
    text=['Controlling Criterion is ', char(CritName), ' Component ',num2str(Comp)];
else
    text=['Controlling Criterion is ', char(CritName)];
end
style='Normal';
WriteText(Format, fid, Fmt2, word, text, style, [1,0])

style='Normal';
text=['Minimum MoS: ', num2str(MoS.MinMoS)];
WriteText(Format, fid, Fmt2, word, text, style, [1,0])


end


