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

        function WitnessOp = calculate_mu_witness(MuTree, Povm, k, s_k)
          % Povm is the whole tree
          WitnessOp = eye(4)*MuTree{k, s_k};
          n_outputs = size(Povm{k, s_k}, 2);
          for j=1:(n_outputs-1)
            % disp(MuTree{k+1, j})
            WitnessOp = WitnessOp - (MuTree{k+1, j})*Povm{k, s_k}{j};
          end
        end
        
        function OptPovmTree = optimize_povm_round(OmTree, PovmTree, StateDims, MuTree, k, KTilde, err1)
          % Optimization function. This updates the POVMs of all states in
          % round k. The SDP is an application of Eq53.
          % Optimization Variables:
          %       (1)  {M_{s'|s}} for all s' in S_{k+1} and s in S_k
          %       (2)  {mu^j_{s_j}} for j<=k  and all s_j
          %
          % Min        Err_II
          % subject to:
          %            Err_II >= 0
          %            Omega - (1 - Err_II)*Eye(d_E) >= 0
          %            for all s in S_k   and s' in S_{k+1}
          %                M_{s'|s} >= 0
          %                Sum_{s'} M_{s'|s} == Eye(4)
          %            for all j<= k , and s in S_j ::
          %                mu^j_{s}*Eye(4) - Sum_{s'} M_{s'|s} * mu^{j+1}_{s'} == W_{j, s}
          %                W_{j, s} = W1_{j,s} + PartialTranspose(W2_{j,s})
          %                W1_{j,s} , W2_{j,s} >= 0
          %                0<= mu^j_{s} <= 1
          %            mu^1_1 <= error_1

          
          n_outputs = size(PovmTree{k,1}, 2);  %
          dE = size(KTilde{1}, 1);
          
          err2 = sdpvar(1,1); % target value to minimize
          constraints = [err2>=0];
          
          % OptPovmTree is like PovmTree except it has SDP variables for the corresponding state
          OptPovmTree = PovmTree;
          for s_k=1:StateDims(k)
              OptPovmTree{k, s_k}{n_outputs} = eye(4);
              for i=1:(n_outputs-1)
                OptPovmTree{k, s_k}{i} = sdpvar(4, 4, 'hermitian', 'complex');
                OptPovmTree{k, s_k}{n_outputs} = OptPovmTree{k, s_k}{n_outputs} - OptPovmTree{k, s_k}{i};
                constraints = [constraints, OptPovmTree{k, s_k}{i}>=0];
              end
              constraints = [constraints, OptPovmTree{k, s_k}{n_outputs}>=0];
          end

          % Create Mu Variables
          OptMuTree = MuTree;
          if k>1
              for s_k=1:StateDims(k)
                  OptMuTree{k, s_k} = sdpvar(1,1);
                  constraints = [constraints, OptMuTree{k, s_k}>=0, OptMuTree{k, s_k}<=1];
              end
              for j=(k-1):-1:2
                  for s_j=1:StateDims(j)
                      OptMuTree{j, s_j} = sdpvar(1,1);
                      constraints = [constraints, OptMuTree{j, s_j}>=0, OptMuTree{j, s_j}<=1];
                  end
              end
          end
          OptMuTree{1, 1} = sdpvar(1,1);
          constraints = [constraints, OptMuTree{1, 1}>=0, OptMuTree{1, 1}<=err1];
          
          % Create Witness Variables and constraints
          for s_k=1:StateDims(k)
              W{k, s_k}{1} = sdpvar(4, 4, 'hermitian', 'complex');
              W{k, s_k}{2} = sdpvar(4, 4, 'hermitian', 'complex');
              constraints = [constraints, W{k, s_k}{1}>=0, W{k, s_k}{2}>=0];
              witness{k, s_k} = CS.calculate_mu_witness(OptMuTree, OptPovmTree, k, s_k);
              constraints = [constraints, witness{k, s_k} == (W{k, s_k}{1}+QI.PT(W{k, s_k}{2},2,[2,2]))]; 
          end
          for j=(k-1):-1:1
              for s_j=1:StateDims(j)
                 W{j, s_j}{1} = sdpvar(4, 4, 'hermitian', 'complex');
                 W{j, s_j}{2} = sdpvar(4, 4, 'hermitian', 'complex');
                 constraints = [constraints, W{j, s_j}{1}>=0, W{j, s_j}{2}>=0];
                 witness{j, s_j} = CS.calculate_mu_witness(OptMuTree, OptPovmTree, j, s_j);
                 constraints = [constraints, witness{j, s_j} == (W{j, s_j}{1}+QI.PT(W{j, s_j}{2},2,[2,2]))];
              end
          end
          
          % Calculate Variable omega
          OptOmTree = OmTree; % variable held here
          for s_k=1:StateDims(k)
            OptOmTree{k,s_k} = CS.calculate_omega_state(OptOmTree, OptPovmTree{k, s_k}, k, KTilde);
          end
          for j=(k-1):-1:1
            for s_j=1:StateDims(j)
              OptOmTree{j, s_j} = CS.calculate_omega_state(OptOmTree, OptPovmTree{j, s_j}, j, KTilde);
            end
          end
          constraints = [constraints, (OptOmTree{1,1} - (1-err2)*eye(dE))>=0];
          
          % Perform SDP
          options = sdpsettings('verbose', 0, 'solver','mosek');
          optimize(constraints, err2, options);
          
          % Placing the optimized values into the output array
          for s_k=1:StateDims(k)
              for i=1:n_outputs
                OptPovmTree{k, s_k}{i} = value(OptPovmTree{k, s_k}{i});
              end
          end
          yalmip('clear'); % clears memory, otherwise process slows down
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