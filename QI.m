classdef QI
    methods(Static)
    function U = rand_unitary(d)
      % d dim Unitary matrix
      H = randn(d) + 1i*randn(d);
      [Q,R] = qr(H);
      R = sign(diag(R));
      R(R==0) = 1;
      U = bsxfun(@times, Q, R.');
    end

    function POVM = rand_povm(d, m)
      % d dim POVM , m elements
      % Output is a cell array:
      % POVM{i} >= 0
      % Sum_{i=1}^{m}  POVM{i} == eye(d)
      k = floor(d/m)+1;
      U = QI.rand_unitary(m*k);
      V = U(:, 1:d);
      for ii = 1:m
        P = zeros(m);
        P(ii, ii) = 1;
        POVM{ii} = ctranspose(V) * kron(P, eye(k)) * V;
      end
    end

    function POVM = optimized_rand_povm(d, m, q)
      % produces a set of random POVMs, such that 
      % max_{sigma in SEP}   Tr( POVM{m} sigma) = q 
      % q is min prob to output SEP(mth output)
      POVM = QI.rand_povm(d, m);
      f = 1;
      if(m==2)
        % empirically, m==2 is the only moment when the below produces f!= 1
        sum_op = eye(d) - POVM{m};
        f = QI.optimize_over_PPT_states_d4(sum_op);
        if(f==1)
          f = 1 - 1e-10;
        end
      end
      eta = q/f;
      for i=1:m
        POVM{i} = POVM{i}*eta;
      end
      POVM{m} = POVM{m} + (1-eta)*eye(d);     
    end

    function rho = rand_state(d)
      % random density matrix of dim d
      H = randn(d) + 1i*randn(d);
      rho = H*H';
      rho = rho/trace(rho);
    end

    function x = PT(rho, sys, dim)
    % partial transpose, dim=[d_1, d_2, ..., d_n]
    % sys = index to be transposed;  sys \in [n]
      n = length(dim);
      d = size(rho);
      perm = [1:2*n];
      perm([n+1-sys,2*n+1-sys]) = perm([2*n+1-sys,n+1-sys]);
      x = reshape(permute(reshape(rho,[dim(end:-1:1),dim(end:-1:1)]),perm),d);
    end

    function y = is_npt(rho)
      % verifies if a 2-qubit state has negative partial transpose
      % 0 = SEP ; 1 = ENT
      assert(size(rho,1)==4, 'Only Implemented for 2 qubit state.');
      y = min(eig(QI.PT(rho,2,[2,2])))<0 ;
    end

    function phi = bell_state(i)
      % Returns a density matrix of a maximally entangled 2 qubit state
      % i = {1,2,3,4} 
      BellMatrix = (1/sqrt(2))*[1, 0, 0, 1; 0, 1, 1, 0; 0, 1, -1, 0; 1, 0, 0, -1];
      phi = BellMatrix(i,:)'*BellMatrix(i,:); 
    end

    function f = optimize_over_PPT_states_d4(operator)
      %     f =  max  trace( operator * sigma)
      % subject to:
      %            sigma >= 0
      %            trace(sigma) == 1 
      %            partial transpose(sigma) >= 0
      %
      % !!  requires operator to be dim 4  
      assert(size(operator,1)==4, 'Only Implemented for 2-qubit state.');
      sigma = sdpvar(4,4,'hermitian','complex');
      constraints = [sigma>=0, QI.PT(sigma,2,[2,2])>=0, trace(sigma)==1];
      options = sdpsettings('verbose', 0, 'solver','mosek');
      optimize(constraints, -trace(operator*sigma), options);
      f = value(trace(operator*sigma));
      yalmip('clear');
    end
    
    function adg = a_dagger(n_levels)
      % Ladder Operator
      % The ladder does not wrap around itself
      % a_dagger |n> = 0
      adg = zeros(n_levels);
      for i=1:(n_levels-1)
        adg(i+1,i) = sqrt(i);
      end
    end

    function H = h_interaction(n_levels)
      % A simple interaction hamiltonian as an example
      % The environment dimension d_E = n_levels
      % System dimension is 4 = 2 qubits
      hE = QI.a_dagger(n_levels);
      sigma_minus = zeros(2);
      sigma_minus(2, 1)=1;
      hAB = kron(eye(2), sigma_minus) + kron(sigma_minus, eye(2));
      H = kron(hE, hAB);
      H = H + H';
    end

    function U = unitary_step(hamiltonian, time)
      % Creates a Unitary that is equivalent of evolving the hamiltonian
      % for a specific amount of time.
      U = expm(-1i*time*hamiltonian);
    end

    function RHO = fin_corr_family(rho_E_0, rho_AB_0, n_rounds, time)
     % RHO is a 1-d cell_array RHO{k} that contains 4x4 matrices
     % rho_E_0 is the initial environment state
     % rho_AB_0 is the ideal target state (must be 4x4)
     % n_rounds is size of output RHO array
     % time is how long each rho_AB_0 interacts with rho_E_t. good value is ~0.1
      dE = size(rho_E_0, 1);
      H = QI.h_interaction(dE);
      U = QI.unitary_step(H, time);
      RHO{n_rounds} = 0; % initialize RHO
      rho_E_t = rho_E_0;
      for t=1:n_rounds
        rho_EAB = kron(rho_E_t, rho_AB_0);
        rho_EAB = U*rho_EAB*U';
        RHO{t} = PartialTrace(rho_EAB, 1, [dE,4]); % from QETLAB
        rho_E_t = PartialTrace(rho_EAB, 2, [dE,4]); % from QETLAB
      end
    end 

    function K = kraus_unitary_step(hamiltonian, time)
      % The user tries to produce the ideal maximally entangled state
      % |PHI+> = (1/sqrt(2))* [1,0,0,1] 
      % But at each round an uncontrolled environment interacts with the
      % ideal state for a given amount of time, and fixed interaction
      % Hamiltonian. This function produces the kraus operator that
      % corresponds to the environment evolution after each round.
      dE = int8(size(hamiltonian, 1)/4);
      Bell = (1/sqrt(2))*[1, 0, 0, 1; 0, 1, 1, 0; 0, 1, -1, 0; 1, 0, 0, -1];
      U = QI.unitary_step(hamiltonian, time);
      U = U * kron(eye(dE), Bell);
      ind2drop = 1:dE*4;
      for i=dE*4:-1:1
        if(mod(ind2drop(i),4)==1)
          ind2drop(i)=[];
        end
      end
      K = U;
      K(:,ind2drop)=[];
      % can check property: K'*K = eye(dE)
    end

    function B = bra_i(i, dim)
      % computational basis bra
      B = zeros(1,dim);
      B(i) = 1;
    end

    function K_tilde = kraus_tilde_interaction(n_levels, time)
      % produces a cell array of size 4, of the main Kraus operator
      % example, tensored computational basis bra of the system qubits:
      % K_tilde{i} = (EYE_{d_E} \otimes <i| ) K
      H = QI.h_interaction(n_levels);
      K = QI.kraus_unitary_step(H, time);
      for i=1:4
        b_i = QI.bra_i(i, 4); 
        K_tilde{i} = kron(eye(n_levels), b_i)*K;
      end
    end
    
    end
end