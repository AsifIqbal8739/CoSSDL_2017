%% Find the atom number with most correlation with the given time series
% Outputs v (correlation value) and d (atom #) Sign is preserved in v
function [v,d] = find_d(Dict,T_series)
    Cor = (corr(Dict,T_series));    Cor(isnan(Cor)) = 0;
    [V,In] = sort(abs(Cor),1,'descend');
    d = In(1); %v = V(1);
    v = Cor(In(1));
end