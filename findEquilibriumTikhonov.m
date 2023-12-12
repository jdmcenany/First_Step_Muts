function [Harvest, Delta] = findEquilibriumTikhonov(prop,tol)
% Adapted from Tikhonov & Monasson 2017, with their CRM (for Fig S2A-B)
% Find ecological equilibrium of a community defined by prop

RANK_STEP = 5000;
R = (prop.capacity(:))';
N = prop.N;
fun = @(h)funToMinimize(h, R);
hessianfcn = @(h,lambda)Hessian(h, R);

if nargin < 2
    tolVal = 1e-11;
else
    tolVal = min(1e-11,tol);
end
%options = optimoptions('fmincon','GradObj','on');
options = optimoptions('fmincon','GradObj','on','algorithm','interior-point',...
    'Hessian','user-supplied','HessFcn', hessianfcn, ...%'DerivativeCheck','on',...
    'TolCon',tolVal,'TolFun',tolVal,'TolX',tolVal,'Display','none');
h0 = 0.9*ones(1,N);
Aeq = [];
beq = [];
lb = zeros(1,N);
ub = Inf(1,N);
nonlcon = [];

% if the number of constraints is very large, optimization is slow
% Thankfully, we know that only a small number of them actually counts!
% So let's only include a small set of species into consideration at first,
% and expand only if necessary. The formula for choosing the initial
% fraction is empirically informed, to work most of the time. But if we
% choose it wrong, we are guaranteed to catch the error, so it's only a
% question of runtime, not accuracy.


A = double(prop.enzymesInSpecies);
b = prop.costs;
% sort by increasing cost
[b, order] = sort(b,'ascend');
A = A(order,:);
maxRank = min(RANK_STEP, length(order));
done = false;
while ~done
    subset = 1:maxRank;
    % when optimizing, use only a subset of constraints...
    Harvest = fmincon(fun,h0,A(subset,:),b(subset),Aeq,beq,lb,ub,nonlcon,options);
    Delta = A*Harvest'-b;
    % ... but verify that ALL the constraints are satisfied in the end
    done = all(Delta<=tolVal*50);
    if ~done
        fprintf('Some ignored constraints are violated; re-trying.\n');
        if maxRank == length(order)
            fprintf(string(max(Delta)));
            Delta = 1000 + Delta;
            break
        end
        maxRank = min(maxRank+RANK_STEP,length(order));
        % increase the subset
    end
end 
Delta(order) = Delta;
% present = abs(Delta)<1e-6;
% 
% DeltaSrt = sort(Delta,'descend');
% semilogy(-(DeltaSrt))
end


function [f,g] = funToMinimize(h, R)
% R is a row vector of resource capacities
% Calculate objective f
f = -R*log(h(:));

if nargout > 1 % gradient required
    g = -R(:)./h(:);
end
end

function H = Hessian(h,R)
    H = diag(R(:)./(h(:).^2));
end
