classdef CS
    methods(Static)
        function OP = omega_product(X, KTilde, PovmM)
          % Define an Omega Product of matrix X and POVM element M:
          % Sum_{ij}   KT_j ' * X * KT_i *(<i|M|j>)
          % where KT = Kraus Tilde operators
          assert(size(PovmM,1)==4, 'Only implemented for 2 qubits.');
          OP = zeros(size(KTilde{1}, 1));
          for i=1:4
            b_i = QI.bra_i(i, 4);
            for j=1:4
              k_j = QI.bra_i(j, 4)';
              OP = OP + (KTilde{j}' * X * KTilde{i})*( b_i*PovmM*k_j);
            end
          end
        end

        function OmState = calculate_omega_state(OmTree, Povm, k, KTilde)
          % Calculates the matrix Omega^j_{s_k}  of round j, state s_k
          % Povm is a 1-d cell array of the game state s_k
          % Assumes that OmTree{k+1} is populated correctly
          OmState = zeros(size(OmTree{k+1,1}, 1));
          n_outputs = size(Povm, 2);
          for i=1:(n_outputs-1)
            OmState = OmState + CS.omega_product(OmTree{k+1,i}, KTilde, Povm{i});
          end
        end

        function OmTree = calculate_omega_tree(PovmTree, StateDims, KTilde)
          % Calculates Omega^k_{s_j}  which are needed to calculate the
          % total operator Omega, from which Error II can be bounded.
          num_mmts = size(StateDims, 2);
          dE = size(KTilde{1}, 1);
          OmTree{num_mmts+1,1} = eye(dE);
          for k=num_mmts:-1:1
            for s_k=1:StateDims(k)
              OmTree{k, s_k} = CS.calculate_omega_state(OmTree, PovmTree{k, s_k}, k, KTilde);
            end
          end
        end

        function OptOmPovm = optimize_omega_state(OmTree, PovmTree, StateDims, MuTree, k, s_k, KTilde)
          % Optimization function. This updates the POVMs after doing an
          % optimization on configuration state {k,s_k}.
          % Minimizes Err_II
          % subject to:
          %            Err_II >= 0
          %            M_i >=0  for all i={1,2,... ,m-1}
          %            W_1>=0 , W_2>=0
          %            Eye(4) - Sum_{i=1}^{m-1} M_i >= 0
          %            mu^k_{s}*Eye(4) - Sum_{i=1}^{m-1} M_i * mu^{k+1}_i >= W_1 + Partial Transpose(W_2)
          %            Omega - (1 - Err_II)*Eye(d_E) >= 0
          n_outputs = size(PovmTree{k,s_k}, 2);
          dE = size(KTilde{1}, 1);
          OptOmTree = OmTree;
          err2 = sdpvar(1,1);
          constraints = [err2>=0];
          for i=1:(n_outputs-1)
            var_povm{i} = sdpvar(4, 4, 'hermitian', 'complex');
            constraints = [constraints, var_povm{i}>=0];
          end
          var_povm{n_outputs} = zeros(4);
          % Calculate Variable omega
          OptOmTree{k,s_k} = CS.calculate_omega_state(OptOmTree, var_povm, k, KTilde);
          for j=(k-1):-1:1
            for s_j=1:StateDims(j)
              OptOmTree{j, s_j} = CS.calculate_omega_state(OptOmTree, PovmTree{j, s_j}, j, KTilde);
            end
          end
          omega = OptOmTree{1,1};
          % Create Witness Variables
          W1 = sdpvar(4, 4, 'hermitian', 'complex');
          W2 = sdpvar(4, 4, 'hermitian', 'complex');
          constraints = [constraints, W1>=0, W2>=0];

          witness = eye(4)*MuTree{k, s_k};
          sum_op = zeros(4);
          for j=1:(n_outputs-1)
              witness = witness - var_povm{j}*MuTree{k+1,j};
              sum_op = sum_op + var_povm{j};
          end
          constraints = [constraints, (eye(4)-sum_op)>=0];
          constraints = [constraints, witness==(W1+QI.PT(W2,2,[2,2]))];
          constraints = [constraints, (omega - (1-err2)*eye(dE))>=0];
          options = sdpsettings('verbose', 0, 'solver','mosek');
          optimize(constraints, err2, options);
          % disp('Error II:')
          % disp(value(err2))
          OptOmPovm = PovmTree;
          OptOmPovm{k, s_k}{n_outputs} = eye(4);
          for i=1:(n_outputs-1)
            OptOmPovm{k, s_k}{i} = value(var_povm{i});
            OptOmPovm{k, s_k}{n_outputs} = OptOmPovm{k, s_k}{n_outputs} - OptOmPovm{k, s_k}{i};
          end
          yalmip('clear');
        end

        function ERROR2 = calculate_error2(OmTree)
          % Calculates the worst case probability that a SEP state is
          % classified as Entangled. This is done assuming that the
          % environment is interacting with the produced quantum states at
          % each round according to the Kraus Operator(s).
          err2 = sdpvar(1,1);
          omega = OmTree{1,1};
          dE = size(omega, 1);
          constraints = [err2>=0, (omega - (1-err2)*eye(dE))>=0];
          options = sdpsettings('verbose', 0);%;, 'solver','sedumi');
          optimize(constraints, err2, options);
          ERROR2 = value(err2);
        end
        
        function MuTree = calculate_mu_tree(PovmTree, StateDims)
          % Creates a mu_tree - necessary to calculate Error Type I
          num_mmts = size(StateDims, 2);
          for k=num_mmts:-1:1
            for s_k=1:StateDims(k)
              if (k==num_mmts)
                % tic
                MuTree{k, s_k} = QI.optimize_over_PPT_states_d4(PovmTree{k,s_k}{1});
                % toc
              else
                % tic
                operator = zeros(4);
                for s_prime=1:StateDims(k+1)
                  operator = operator + PovmTree{k,s_k}{s_prime}*MuTree{k+1, s_prime};
                end
                MuTree{k, s_k} = QI.optimize_over_PPT_states_d4(operator);
                % toc
              end
            end
          end
          MuTree{num_mmts+1,1} = 1;
        end
        
        function PovmTree = init_povm_tree(StateDims, Err1)
          % Initializes a set of configuration states according to
          % StateDims. At each point in the tree indexed as {k, s_k} for
          % round k, state s_k : there is a set of POVMs. Furthermore the
          % last element m in a state is the "Output to SEP" POVM. The
          % initialization optimizes the whole tree such that the Type I
          % error is Err1. I.e. mu{1,1}=err1
          assert(StateDims(1)==1) % only 1 initial state
          dAB = 4;     % dAB = 4 default
          num_mmts = size(StateDims, 2);
          min_prob_sep = Err1^(1/num_mmts);
          % Creates a povm_tree
          for k=1:num_mmts
            for s_k=1:StateDims(k)
              if (k~=num_mmts)      
                PovmTree{k, s_k} = QI.optimized_rand_povm(dAB, StateDims(k+1)+1, min_prob_sep);
              else
                PovmTree{k, s_k} = QI.optimized_rand_povm(dAB, 2, min_prob_sep);
              end
            end
          end
          % povm_tree{k,s}{s'} is a 4x4 matrix povm element
          % sum_s'  povm_tree{k,s}{s'} = eye(4)
        end
        
        function NewMuTree = update_mu_tree(MuTree, PovmTree, StateDims, k, s_k)
          % If only the POVMs of state s_k in round k were changed, then
          % this is an efficient way of updating the mu tree without
          % recalculating everyting
          % see update_mu_tree_round if all states in round k were updated
          NewMuTree = MuTree;
          num_mmts = size(StateDims, 2);
          if(k==num_mmts)
            NewMuTree{k, s_k} = QI.optimize_over_PPT_states_d4(PovmTree{k,s_k}{1});
          else
            operator = zeros(4);
            for s_prime=1:StateDims(k+1)
              operator = operator + PovmTree{k,s_k}{s_prime}*MuTree{k+1, s_prime};
            end
            NewMuTree{k, s_k} = QI.optimize_over_PPT_states_d4(operator);
          end
          for j=(k-1):-1:1
            for s_j=1:StateDims(j)
              operator = zeros(4);
              for s_prime=1:StateDims(j+1)
                operator = operator + PovmTree{j,s_j}{s_prime}*NewMuTree{j+1, s_prime};
              end
              NewMuTree{j, s_j} = QI.optimize_over_PPT_states_d4(operator);
            end
          end
        end
        
        function NewMuTree = update_mu_tree_round(MuTree, PovmTree, StateDims, k)
          %if some or all POVMs of configuration states of round k were
          %changed, then this function updates the tree without needing to
          %recalculate everything
          NewMuTree = MuTree;
          num_mmts = size(StateDims, 2);
          if(k==num_mmts)
            NewMuTree = CS.calculate_mu_tree(PovmTree, StateDims);
          else
            for j=k:-1:1
              for s_j=1:StateDims(j)
                operator = zeros(4);
                for s_prime=1:StateDims(j+1)
                  operator = operator + PovmTree{j,s_j}{s_prime}*MuTree{j+1, s_prime};
                end
                NewMuTree{j, s_j} = QI.optimize_over_PPT_states_d4(operator);
              end
            end
          end
        end
        
    end
end