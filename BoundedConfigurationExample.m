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


n_levels = 30; % environment dimension
interaction_time = 0.1; % How long each |Phi_+> interacts with the environment
% setting the interaction time to 0 is the idealized scenario.
KT = QI.kraus_tilde_interaction(n_levels, interaction_time); % Kraus ops.

num_measurements = 10; % Number of game rounds
bounded_size = 5; % Size of |S_k| excluding the SEP output.
%i.e.  bounded_size + 1 = |S_k| 
state_dims = ones(1, bounded_size)*num_measurements;
state_dims(1)=1;

num_epochs = 5;
% How many times each set of POVMs is optimized.
ERR{num_epochs} = [[,]]; % plot points stored here


% This file generates a plot with Type I Error (e1) on the x-axis:  
for error_1=0:0.5:1
    % for each error_1, create a graph (see Fig 10 in paper) that contains
    % POVM elements for each game state. After that, calculate mu (eq 8)
    % and Omega (eq 12) recursively. 
    povm_tree = CS.init_povm_tree(state_dims, error_1);
    mu_tree = CS.calculate_mu_tree(povm_tree, state_dims);    
    omega_tree = CS.calculate_omega_tree(povm_tree, state_dims, KT);
    for epoch=1:num_epochs
      for k=1:size(state_dims,2)
        for s_k=1:state_dims(k)
          % for each state s_k of game round k, optimize the POVMs
          povm_tree = CS.optimize_omega_state(omega_tree, povm_tree, state_dims, mu_tree, k, s_k, KT);
          % after each POVM optimization, the recursive Omega^(k)_{s_k} need to be updated
          omega_tree = CS.calculate_omega_tree(povm_tree, state_dims, KT);
        end
      end
      % finally, after 1 epoch, the mu's can be updated
      mu_tree = CS.calculate_mu_tree(povm_tree, state_dims);
      new_errors=[mu_tree{1,1}, CS.calculate_error2(omega_tree)];
      disp(new_errors);
      ERR{epoch}=[ERR{epoch};new_errors];
    end
end

% at the end, the plot points are stored in ERR{num_epochs}

plot(ERR{num_epochs}(:,1),ERR{num_epochs}(:,2),'.','MarkerSize',15)
xlabel('e_I', 'Interpreter', 'tex', 'FontSize', 20) 
ylabel('e_{II}', 'Interpreter', 'tex', 'FontSize', 20) 
