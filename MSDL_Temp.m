%% Mult-Subject Dictionary learning along with Common Dictionary for Temporal Case
% [D0,Di,X0,Xi] = argmin 0.5||Yi - D0X0 - DiXi||_F^2 + eta*||Di*A_/i||F^2
% In this script, we learn X0 -> Xi -> D0 -> Di

function [Dict_0,Dict,X_0,X] = MSDL_Temp(St)
%% Some Data Setup
[n,N] = size(St.Yn);    Ns = N/St.nSub;
Data = cell(1,St.nSub);
for i = 1:St.nSub
    Data{i} = St.Yn(:,(i-1)*Ns+1:i*Ns);
end

%% Dictionary Initialization
X = cell(1,St.nSub);                                  % Sub Specific Sparse Codes
X_0 = zeros(St.Size_D0,Ns);
Dict = cell(1,St.nSub);
for i = 1:St.nSub
    if St.DiniType == 1
        if i == 1
            Dict_0 = normc(St.Yn(:,randperm(N,St.Size_D0)));
        end
        Dict{i} = normc(St.Yn(:,randperm(N,St.Size_Ds)));
    else
        if i == 1
            Dict_0 = normc(randn(n,St.Size_D0));
        end
        Dict{i} = normc(randn(n,St.Size_Ds));
    end
    X{i} = zeros(St.Size_Ds,Ns);
end

%% Dict Learning
for iter = 1:St.nIter
	for i=1:2
		[X_0,X] = SparseCode(Data,Dict_0,Dict,X_0,X,St.Spar,St.AlgoType);
	end
    if St.verbose == 1; displayError(Data,Dict_0,Dict,X_0,X,iter);  end
    
    [Dict_0,Dict] = DictUpdate2(Data,Dict_0,Dict,X_0,X,St);
    if St.verbose == 1; displayError(Data,Dict_0,Dict,X_0,X,iter);  fprintf('\n');  end
    
end
% fprintf('............END..........\n');   
end

%% Sparsity Implementation for MSDL_Temporal Scene
% Old Sparse codes are needed for this step
function [X_0,X] = SparseCode(Data,Dict_0,Dict,X_0,X,Spar,AlgoType)    
    M = size(Dict,2);
    
    if AlgoType == 1
        X_0 = SparseCommon;
        X   = SparseIndividual;
    elseif AlgoType == 2
        X   = SparseIndividual;
        X_0 = SparseCommon;
    end        
        
    
    % Inner Functions
    function X_0_ = SparseCommon
        Et = 0;
        for i = 1:M
            Et = Et + (Data{i} - Dict{i}*X{i}); 
        end
        Et = Et./M;
        X_0_ = omp(Dict_0,Et,Dict_0'*Dict_0,Spar(1));
    end

    function X_ = SparseIndividual
        for i = 1:M
            Ei = Data{i} - Dict_0*X_0;
            X_{i} = omp(Dict{i},Ei,Dict{i}'*Dict{i},Spar(2));
        end
    end
end

%% Dictionary Update Step Gradient Descent
function [Dict_0,Dict] = DictUpdate(Data,Dict_0,Dict,X_0,X,St)
    M = St.nSub;
    if St.AlgoType == 1
        DictCommon;
        DictIndividual;
    elseif St.AlgoType == 2
        DictIndividual;
        DictCommon;
    end

    % Inner Functions
    function DictCommon
        W = 0;  A0 = [];
        for i = 1:M    
            A0 = [A0,Dict{i}];
            W = W + (Data{i} - Dict{i} * X{i});
        end
        DD = Dict_0;    W = W./M;   
        A00 = A0*A0';   rho = 1/norm(full(X_0)*full(X_0)');
        for k = 1:St.Kmax
            tt = rho*((DD*X_0-W)*X_0' + St.eta*(A00)*DD);
            DD = normc(DD - tt);    %- eye(size(DD,1))
            if norm(tt,'fro') < St.eps;         break;   end
        end
        Dict_0 = DD;
    end

    function DictIndividual
        for i = 1:M                
            Ai = Dict_0;
            for h = 1:M   
                if h == i; continue; end 
                Ai = [Ai,Dict{h}];   
            end;
            Gi = Data{i} - Dict_0*X_0;
            DD = Dict{i};   XX = X{i};
            Aii = Ai*Ai';   rho = 1/norm(full(XX)*full(XX)');
            for k = 1:St.Kmax
                tt = rho*((DD*XX-Gi)*XX' + St.eta*(Aii)*DD);    % - eye(size(DD,1))
                DD = normc(DD - tt);
                if norm(tt,'fro') < St.eps;    break;   end
            end
            Dict{i} = (DD);
        end        
    end
end

%% Approximation Error to display
function displayError(Data,Dict_0,Dict,X_0,X,iter)
    for b = 1:size(Dict,2)
       EE(b) = (norm(Data{b} - Dict_0*X_0 - Dict{b}*X{b},'fro')); 
    end
    display([sprintf('Iter:%2d, Total Error: ',iter),sprintf('%0.2f,\t',EE)]) %
end

%% Dictionary Update Step Using ADMM
function [Dict_0,Dict] = DictUpdate2(Data,Dict_0,Dict,X_0,X,St)
    M = St.nSub;    max_mu = St.max_mu;     scale_mu = St.scale_mu;   ini_mu = St.ini_mu;
    Z0 = zeros(size(Dict_0));   mu = ini_mu;    
    Zi = zeros(size(Dict{1}));  W0 = Z0;    Wi = Zi; 
    if St.AlgoType == 1
        DictCommon;
        DictIndividual;
    elseif St.AlgoType == 2
        DictIndividual;
        DictCommon;
    end

    % Inner Functions
    function DictCommon
        T = 0;  A0 = [];
        Z0 = zeros(size(Dict_0));   mu = ini_mu; W0 = Z0;
        for i = 1:M    
            A0 = [A0,Dict{i}];
            T = T + (Data{i} - Dict{i} * X{i});
        end
        DD = Dict_0;    T = T./M;   
        A00 = A0*A0';   %rho = 1/norm(full(X_0)*full(X_0)');
        for k = 1:St.Kmax
%             while norm((DD - Z0),'fro') > 0.01
                DD = (mu*Z0 + T*X_0' - W0)/(X_0*X_0' + mu*eye(size(DD,2)));  DD = normc(DD);
                Z0 = (mu*eye(size(DD,1))+ 2*St.eta*A00)\(W0 + mu*DD); Z0 = normc(Z0);
%             end
            W0 = W0 + mu*(DD-Z0);
            mu = min(scale_mu*mu,max_mu);
%             display(norm((DD - Z0),'fro'))
            if norm((DD - Z0),'fro') < St.eps;  break;  end
        end
        Dict_0 = DD;
    end

    function DictIndividual
        for i = 1:M        
            Zi = zeros(size(Dict{i}));   mu = ini_mu; Wi = Zi;
            Ai = Dict_0;
            for h = 1:M   
                if h == i; continue; end 
                Ai = [Ai,Dict{h}];   
            end;
            Gi = Data{i} - Dict_0*X_0;
            DD = Dict{i};   XX = X{i};
            Aii = Ai*Ai';   % rho = 1/norm(full(XX)*full(XX)');
            for k = 1:St.Kmax
%                 while norm((DD - Zi),'fro') > 0.01
                    DD = (mu*Zi + Gi*XX' - Wi)/(XX*XX' + mu*eye(size(DD,2)));  DD = normc(DD);
                    Zi = (mu*eye(size(DD,1))+ 2*St.eta*Aii)\(Wi + mu*DD); Zi = normc(Zi);
%                 end
                Wi = Wi + mu*(DD-Zi);
                mu = min(scale_mu*mu,max_mu);
            
                if norm((DD - Zi),'fro') < St.eps;    break;         end
            end
            Dict{i} = (DD);
        end        
    end
end






