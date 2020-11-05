%-------------------------------------------------------------------------%
%       Example for bounded configuration space entanglement detection
%
%       This script uses Yalmip (with the solver MOSEK)
%       And the auxiliary function files QI.m and CS.m
%-------------------------------------------------------------------------%

% In this example, the user wants to create a maximally entangled state of
% 2 qubits in each round. (This is encoded in QI.kraus_unitary_step).
% Before leaving the player's lab, the state interacts with an environment
% briefly.

% RUNTIME: On a normal laptop - Around 480 sec. 
% For a quicker demonstration, consider lowering the num_measurements.

% ----------- Environment Parameters ---------
% --------------------------------------------
n_levels = 10; % environment dimension
interaction_time = 0.05; % How long each |Phi_+> interacts with the environment
% setting the interaction time to 0 is the idealized scenario.
KT = QI.kraus_tilde_interaction(n_levels, interaction_time); % Kraus ops.

% ------------ Game Settings -----------
% --------------------------------------
num_measurements = 6; % Number of game rounds
bounded_size = 3; % (bounded_size + 1) == |S_k| == Number of outcomes per measurement

% ---- misc. ----
state_dims = ones(1, num_measurements)*bounded_size;
state_dims(1)=1;
% This is how many times each round is optimized
num_epochs = 3;
ERR{num_epochs} = [[,]]; % plot points stored here
disp(['     err1   ','   err2'])

% Main loop for Plot
% (x,y) = (e_I, e_II) 
for error_1 = 0.05:0.1:0.95 % 0.01:0.06:1
    % for each error_1, create a graph (see Fig 10 in paper) that contains
    % POVM elements for each game state. After that, calculate mu (eq 8)
    % and Omega (eq 12) recursively. 
    povm_tree = CS.init_povm_tree(state_dims, error_1);
    mu_tree = CS.calculate_mu_tree(povm_tree, state_dims);    
    omega_tree = CS.calculate_omega_tree(povm_tree, state_dims, KT);
    for epoch=1:num_epochs
        for k=size(state_dims,2):-1:1
          % for each game round k, optimize the POVMs
          % after each POVM round optimization, the recursive Omega^(k)_{s_k} and mu^(k)_{s_k} need to be updated
          povm_tree = CS.optimize_povm_round(omega_tree, povm_tree, state_dims, mu_tree, k, KT, error_1);          
          omega_tree = CS.calculate_omega_tree(povm_tree, state_dims, KT);
          mu_tree = CS.update_mu_tree_round(mu_tree, povm_tree, state_dims, k);
        end
        new_errors=[mu_tree{1,1}, CS.calculate_error2(omega_tree)];
        ERR{epoch}=[ERR{epoch};new_errors];
    end
    disp(new_errors)
end

% At the end, the plot points are stored in ERR{num_epochs}
% Because this is a see-saw method by rounds, the optimization is
% stochastic. Sampling multiple times will yield different results.
plot(ERR{num_epochs}(:,1),ERR{num_epochs}(:,2),'.','MarkerSize',15)
xlabel('e_I', 'Interpreter', 'tex', 'FontSize', 20) 
ylabel('e_{II}', 'Interpreter', 'tex', 'FontSize', 20) 
