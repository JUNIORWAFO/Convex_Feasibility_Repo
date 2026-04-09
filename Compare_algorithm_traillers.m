clear all
clc
addpath("./Utils")
L = 2;

p = gcp('nocreate');
if isempty(p)
    % There is no parallel pool
    poolsize = 0;
    p = parpool('Processes');
else
    % There is a parallel pool of <p.NumWorkers> workers
    poolsize = p.NumWorkers;
end

Nmax = 3;   % max number of trailers
% Results matrices:
Times_tzanis = zeros(Nmax,1);         % stores running times
Times_feas = zeros(Nmax,1);         % stores running times
Volumes_feas = zeros(Nmax,1);       % stores computed volumes
feas_check = zeros(Nmax,1); %Store 
%%%%%%%%%%%%%%%%%%%%%%%% Model parameters %%%%%%%%%%%%%%%%%%%%%%%%
kd = 4600;  % stiffness
ks = 4500;  % damper coefficient
m0 = 500;   % weight of truck
m = 1000;   % weight of trailer
T = 0.4;    % sampling rate

for N = 1:4      % number of trailers
    disp(N);
    
    % Continuous-time system - x = [d1 .. dN v0 v1 .. vN] :
    tmpA12 = [eye(N) zeros(N,1)]+[zeros(N,1) -eye(N)];
    tmpA1 = [zeros(N,N) tmpA12];
    tmpA2 = [-ks/m0 zeros(1,N-1) -kd/m0 kd/m0 zeros(1,N-1)];
    tmpA31 = ks/m*([eye(N-1) zeros(N-1,1)]+[zeros(N-1,1) -eye(N-1)]);
    tmpA32 = kd/m*([eye(N-1) zeros(N-1,1) zeros(N-1,1)]+[zeros(N-1,1) -2*eye(N-1) zeros(N-1,1)]+[zeros(N-1,1) zeros(N-1,1) eye(N-1)]);
    tmpA3 = [tmpA31 tmpA32];
    tmpA4 = [zeros(1,N-1) ks/m zeros(1,N-1) kd/m -kd/m];
    tmpA = [tmpA1; tmpA2; tmpA3; tmpA4];
    tmpB = [zeros(N,1); 1; zeros(N,1)];
    
    % Discrete-time system:
    % Use forward Euler methnod to obtain discrete-time system, with sampling time T:
    Ao = T*tmpA + eye(2*N+1);
    Bo = T*tmpB;
    
    % State constraints:
    d = 0.5;    % spring elongation constraint
    vMax = 10;  % max velocity
    vMin = 0;   % -min velocity
    
    Go =[eye(N) zeros(N,N+1); -eye(N) zeros(N,N+1); zeros(N+1,N) eye(N+1); zeros(N+1,N) -eye(N+1)];
    Fo = [d*ones(N,1); d*ones(N,1); vMax*ones(N+1,1); vMin*ones(N+1,1)];
    D = Polyhedron('H',[Go Fo]);
    h = min(2*N+1,2);
    C = ones(h,size(Go,2));
    C(1:h,1:h) = eye(h);
    Y = C*D;
    Y.computeHRep();
    Y.computeVRep();
    % Final matrices:
    A = Ao; B = Bo; Gx = Go; Fx = Fo;
    U = Polyhedron('V', [0], 'R', [-1;1]);

    n_iter = 5*(N==3) + 0*(N==4);
    n = 2*N+1;
    n_max = 30;
    n_min = h+1;
    n_s = 15;
    iter = 2*N+2;
    %method = 1 *(N<=2) + 2*(N>2);
    method = 0;
    strict = 2;
    tic
    CI_feas(N) = compute_CI_feas_Obs(A,B,C, D,U,Y,n_max,n_min,n_s,23,2,2, 0, 1,'fmincon');
    %CI_feas(N) = backward_compute_3(CI_feas(N),D,U,A,B,n_iter);
    Times_feas(N) = toc;
    CI_feas(N).minVRep();
    disp(Times_feas(N));
    % tic
    %CI_feas(N) = compute_CI_feas(A,B, D,U,n_max,n_min,n_s, n_iter,3,1, method, strict, 'ipopt');
    % CI_feas(N) = backward_compute_3(CI_feas(N),D,U,A,B,n_iter);
    % Times_feas(N) = Times_feas(N)+toc;
    % if N>=3
    %     s = norm(CI_feas(N).V,Inf);
    %     test_vert = round(CI_feas(N).V,N+1);
    %     test_ci = Polyhedron('V', test_vert);
    % end 
    try
        Volumes_feas(N) = CI_feas(N).volume;
        disp(Volumes_feas(N))
    catch
        try
            disp("First volume calculation not working");
            Volumes_feas(N) = test_ci.volume;
            disp(Volumes_feas(N))
        catch
            disp("Volume Calculation not working");
        end
    end
    feas_check(N) = check_invariance(CI_feas(N),A,B,D,U);
    
    save(sprintf('Poly_%d.mat',N), 'CI_feas')
    clear CI_feas
end
% for N =1:Nmax
%     try
%         Volumes_feas(N) = CI_feas(N).volume;
%         disp(Volumes_feas(N))
%     catch
%         try
%             disp("First volume calculation not working");
%             Volumes_feas(N) = test_ci.volume;
%             disp(Volumes_feas(N))
%         catch
%             disp("Volume Calculation not working");
%         end
%     end
%     feas_check(N) = check_invariance(CI_feas(N),A,B,D,U);       
% end