%-------------------------------------------------------------------------%
%       Example for 1-shot entanglement detection
%
%       This script uses CVX (with the solver MOSEK)
%       QETLAB PartialTranspose [http://www.qetlab.com/PartialTranspose]
%       And the auxiliary function file paulis.m
%-------------------------------------------------------------------------%

% Runtime on laptop ~180 sec

% This file generates a plot with Type I Error (e1) on the x-axis:
e1 = 0.01 : 0.01 : 1 ;

% Specify the two entangled states that will be certified.

rho1 = [        0,    0,        0,      0;
                0,  0.5,     -0.5,      0;
                0, -0.5,      0.5,      0;
                0,    0,        0,      0];
 
rho2 = [      0.5,    0, 0.353553, 0.353553;
                0,    0,        0,        0;
         0.353553,    0,     0.25,     0.25;
         0.353553,    0,     0.25,     0.25];
     

n=size(rho1,1);
% Only the first level of the hierarchy is implemented, so results are
% tight only for n=4, otherwise they are an upper bound.

num_pts = size(e1,2);
% Type II Error (e2) optimized for different strategies:
% M_1 : Global Measurements
e2_M1 = zeros(1, num_pts);
% M_2 : LOCC Measurements
e2_M2 = zeros(1, num_pts);
% M_3 : Fixed Measurements
e2_M3 = zeros(1, num_pts);

% Specify the measurements
S1 = paulis(1); % Used for LOCC Measurements
S2 = paulis(2); % Used for Fixed Measurements
  
 for i=1:num_pts

     % Global Measurements: M_1
     cvx_begin sdp
     cvx_solver mosek
                variable x(n,n) hermitian semidefinite 
                variable V1(n,n) hermitian semidefinite
                variable V2(n,n) hermitian semidefinite
                variable e2
                minimize( e2 ); 
                subject to
                    eye(n) - x >= 0
                    e2 >= 1-trace(rho1*x)
                    e2 >= 1-trace(rho2*x)     
                    e1(i)*eye(n)-x == V1 + PartialTranspose(V2,2)
      cvx_end
      e2_M1(i) = cvx_optval

      % LOCC Measurements: M_2
      nin=3;
      nout=2;
      cvx_begin sdp
      cvx_solver mosek 
                variable P(nin,nout,nin,nout,2) 
                variable V1(n,n) hermitian semidefinite
                variable V2(n,n) hermitian semidefinite
                variable e2
                expression P2(nin,nout,nin,nout) 
                expression P3(nin,nout,1,nout)
                expression QQ(n,n) 
                minimize( e2 );
                subject to
                    %LOCC decision constraints
                    for v1=1:nout
                        for v2=1:nout
                            sum(reshape(P(:,v1,:,v2,:),[nin^2*2,1])) == 1 
                        end
                    end
                    reshape(P,[nin^2*nout^2*2,1]) >= 0
                    P2 = sum(P,5);
                    P3 = sum(P2,3);
                    for v3=1:nin
                        for v1=1:nout
                            for v2=1:nout
                                P3(v3,v1,1,v2) == P3(v3,1,1,1)
                                for v4=1:nin
                                    P2(v3,v1,v4,v2) == P2(v3,v1,v4,1)
                                end
                            end
                        end
                    end

                    % separability constraints
                    QQ = zeros(n,n);
                    for v1=1:nin
                        for v3=1:nin
                            for v2=1:nout
                                for v4=1:nout
                                    QQ = QQ + P(v1,v2,v3,v4,1)*kron(cell2mat(S1(v1,v2)),cell2mat(S1(v3,v4)));
                                end
                            end
                        end
                    end
                    e1(i)*eye(n) - QQ == V1 + PartialTranspose(V2,2)

                    % introducing optimisation variable
                    e2 >= 1-trace(QQ*rho1)
                    e2 >= 1-trace(QQ*rho2)
      cvx_end
      e2_M2(i) = cvx_optval

      % Fixed Measurements: M_3
      nin=9;
      nout=4;
      cvx_begin sdp
      cvx_solver mosek
                variable P(nin,nout,2) 
                variable V1(n,n) hermitian semidefinite 
                variable V2(n,n) hermitian semidefinite
                variable e2
                expression P2(nin,nout) 
                expression QQ(n,n) 
                minimize( e2 );
                subject to
                    %LOCC decision constraints
                    for v1=1:nout
                        sum(reshape(P(:,v1,:),[nin*2,1])) == 1 
                    end
                    reshape(P,[nin*nout*2,1]) >= 0
                    P2 = sum(P,3);
                    for v1=1:nin
                        for v3=1:nout
                            P2(v1,v3) == P2(v1,1)
                        end
                    end

                    % separability constraints
                    QQ = zeros(n,n);
                    for v3=1:nin
                        for v4=1:nout
                            QQ = QQ + P(v3,v4,1)*cell2mat(S2(v3,v4));
                        end
                    end
                    e1(i)*eye(n) - QQ == V1 + PartialTranspose(V2,2)

                    % introducing optimisation variable
                    e2 >= 1-trace(QQ*rho1)
                    e2 >= 1-trace(QQ*rho2)
     cvx_end
     e2_M3(i) = cvx_optval
 end


plot(e1(:), e2_M1(:),'.','MarkerSize',10)
hold('on');
plot(e1(:), e2_M2(:),'.','MarkerSize',10)
plot(e1(:), e2_M3(:),'.','MarkerSize',10)
xlabel('e_I', 'Interpreter', 'tex', 'FontSize', 20)
ylabel('e_{II}', 'Interpreter', 'tex', 'FontSize', 20)
legend({'Global','LOCC','Fixed'},'Location','northeast')
hold('off');