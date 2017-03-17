% Generating multi-subject dataset using SIMTB
rng(88);

%% Subject components and other setup
nSub = 6;
IDCommon = [16,23,22];
Comp{1} = [17];      % one SM per subject (these are summed later)
Comp{2} = [19];
Comp{3} = [20];
Comp{4} = [24];
Comp{5} = [7];
Comp{6} = [29,30];     % 1-> Common, others subject specific

sP.nT = 150;
sP.nV = 100;

SubData = cell(1,nSub);

%% SM and TC generation
Mask = simtb_createmask(sP);    Mask = Mask(:);
% Common SMs 
SMCommon = zeros(length(IDCommon),sP.nV^2);
for i = 1:length(IDCommon)
    temp = simtb_generateSM(IDCommon(i),sP.nV);
%     imagesc(temp); title(sprintf('nC:%d',IDCommon(i))); pause
    SMCommon(i,:) = Mask.*temp(:);
end

% Subject Specific SMs
SMSpec = SubData;
for s = 1:nSub
    for i = 1:length(Comp{s})
        temp = simtb_generateSM(Comp{s}(i),sP.nV);
%         imagesc(temp); title(sprintf('SubID:%d, nC:%d',s,Comp{s}(i))); pause
        SMSpec{s}(i,:) = Mask.*temp(:);
    end
    SMSpec{s} = sum(SMSpec{s},1);
end

% TC generation
nTC = 9;       % Count nC manually
TC = zeros(sP.nT,nTC);
for i = 1:nTC
    eTC = zeros(1,sP.nT);
    eTC(randperm(sP.nT-10,20)) = 1;
    TC(:,i) = simtb_TCsource(eTC,2,1);
%     plot(TC(:,i)); pause
end

%% Subject Specific Data Matrix Generation
IDC = length(IDCommon);
TCCommon = TC(:,1:IDC);
DataCommon = TCCommon*SMCommon;     % Common Part

% Subject Data Matrix Scene
for i = 1:nSub
    SubData{i} = DataCommon + TC(:,IDC+i)*SMSpec{i};
%     mesh(SubData{i});pause   
end

%% Saving imp stuff in a structure
Sim.Data = cell2mat(SubData);
Sim.TC = TC;
Sim.SMCommon = SMCommon;
Sim.SMSpec = SMSpec;
Sim.SubData = SubData;
Sim.nT = sP.nT;
Sim.nV = sP.nV;
Sim.nSub = nSub;
Sim.Mask = Mask';

% save('SimTB_Data.mat','Sim');