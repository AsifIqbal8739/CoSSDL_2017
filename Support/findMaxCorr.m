%% Support function for MSDL correlation checks
% TC - GrTr time courses
% SM - GrTr maps
% DD/XX - subject dicts/codes concatenated with common one at the end
function [TCC,SMC] = findMaxCorr(TC,SM,DD,XX)
nComp = size(TC,2);
nS = size(DD,2);  % Last one is common pair
[TCC,SMC] = deal(zeros(nComp,nS));
for iC = 1:nComp
    for iS = 1:nS
        [C,~] = find_d(DD{iS},TC(:,iC));
        TCC(iC,iS) = abs(C);
        [C,~] = find_d(XX{iS}',SM(iC,:)');
        SMC(iC,iS) = abs(C);
    end
end


end