function [A, B, C, R] = cost_function_param(N, X0)
    dt = 60 / N;
    
    A11 = eye(N) - diag(ones(N-1,1),-1);
    A12 = -dt*eye(N);
    A1 = [A11, A12];
    A2 = [zeros(1,N-1), 1, zeros(1,N)];
    A = [A1; A2];
    
    C = [0.1*dt*diag(ones(N-1,1),-1), zeros(N); zeros(1,2*N)];
    
    B1 = [X0 - 0.1*sqrt(X0)*dt; zeros(N-1,1)];
    B2 = 5;
    
    B = [B1; B2];%; B4];
    
    
    R = blkdiag(zeros(N),eye(N));
end