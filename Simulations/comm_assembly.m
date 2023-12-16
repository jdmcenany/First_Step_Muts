function out = comm_assembly(task_id, N_runs, input1, input2, R0, R, input_mode, sigma_R_sq,id_out)
% Randomly assembles a community without mutation
% If input_mode = 0, input1 = S/R and input2 = epsilon
% If input_mode = 1, input1 = <S*/R> and input2 = <S*/S>

function I = I_func(x)
    I = (1+x.^2).*erfc(x/sqrt(2))/2 - x.*exp(-x.^2/2)/sqrt(2*pi);
end

function dI = dI_func(x)
    dI = x.*erfc(x/sqrt(2)) - exp(-x.^2/2)*sqrt(2/pi);
end

rng(task_id)
organisms = cell(N_runs,1);
harvests = cell(N_runs,1);
deltas = cell(N_runs,1);

if input_mode == 1
    alpha_star = input1;
    phi = input2;
    lmbda = sqrt(2)*erfinv(1-2*phi);
    alpha = alpha_star/phi;
    V_tot_root = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
    if V_tot_root < 0
        fprintf('Invalid params: decrease alpha_star or phi.\n');
        return
    end
    epsilon = V_tot_root.^2*(1 - alpha.*I_func(lmbda)) - sigma_R_sq*R0*(1-R0/R)*(1-alpha_star).^2;
    if epsilon < 0
        fprintf('Invalid params: decrease sigma_R_sq, alpha_star or phi.\n');
        return
    end
    epsilon = sqrt(epsilon);
elseif input_mode == 0
    alpha = input1;
    epsilon = input2;
else
    fprintf('Invalid input mode.\n');
    return
end

for i = 1:N_runs
    seed = round(rand*100000000);
    p_val = R0/R;
    R2_val = sigma_R_sq*R;
    
    p.N = R;
    p.R2 = R2_val;
    p.sparsity = p_val;
    rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
    p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
    while sum(p.capacity < 0) > 0
        fprintf('Negative resource supply generated, retrying.\n');
        p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
    end
    p.alpha = alpha;
    p.seed = round(1/epsilon+seed);
    p.epsilon = epsilon/R0;
    prop = generateOrganisms(p);
    prop.capacity = p.capacity;
    true_deltas = {};
    
    [harvest, delta] = findEquilibrium(prop);
    
    organisms{i} = prop;
    harvests{i} = harvest;
    deltas{i} = delta;

    % Check for optimization failure
    if delta(1) > 10 
        out = -3
        return
    end
end

if input_mode == 0
    fname = "Results/comm_alpha_" + string(alpha) + "_epsilon_" + string(epsilon) + "_R0_" + string(R0) + "_R_" + string(R) + "_sigmaR2_" + string(sigma_R_sq) + "_" + string(task_id - id_out) + ".mat";
elseif input_mode == 1
    fname = "Results/comm_astar_" + string(alpha_star) + "_phi_" + string(phi) + "_R0_" + string(R0) + "_R_" + string(R) + "_sigmaR2_" + string(sigma_R_sq) + "_" + string(task_id - id_out) + ".mat";
end

save(fname,'harvests','deltas','organisms','-v7.3')
out = 0;
end
