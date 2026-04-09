function RCI=compute_CI_feas_Obs(A,B,C, X,U,Y,N_max,N_min,n_s,n_iter,n_check,L_check, method, strict,solver)
    rng('shuffle');
    sA = size(A);
    n =  size(A,1);
    m = size(C,1);
    sB = size(B);
    sC = size(C); 
    n_x = size(X.V,1);
    V_X = X.V';
    n_y = size(Y.V,1);
    V_Y = Y.V';
    R = polytope.generateRandom('Dimension',m,'NrConstraints',3*m,...
            'IsDegenerate',false,'IsBounded',true);
    D = R.A;
    R_1 = polytope.generateRandom('Dimension',n,'NrConstraints',3*n,...
            'IsDegenerate',false,'IsBounded',true);
    D_1 = R_1.A;
    test = 0;
    V_null = null(C);
    H_poly = [];
    Z_poly = [];
    c = 0;
    if strcmp(solver, 'fmincon') 
        ops = sdpsettings('solver','fmincon','FinDiffType', 'central','MaxFunEvals',5000, ...
           'verbose',0,'Diagnostics','off','InitBarrierParam', 0.5, 'SubproblemAlgorithm', 'cg');
    elseif strcmp(solver, 'ipopt') 
        ops = sdpsettings('solver','ipopt','mu_strategy', 'monotone', 'verbose',0);
    end
    for N = N_min:N_max
        disp(N);
        x = sdpvar(sA(1),N+1, 'full');
        z = sdpvar(sC(1),N+1, 'full');
        u = sdpvar(sB(2), N, 'full');
        r =sdpvar(3*n,1);
        r_2 =sdpvar(3*m,1);
        Cond = [];
        if strict == 1
            d = sdpvar(2*m,1);
            a = sdpvar(2*m,1);
            lambda_feas = sdpvar(N,2*m, 'full');
            for i = 1:2*m
                h = floor((i-1)/m);
                v = zeros(m,1);
                s = mod(i-1, m) + 1;
                v(s) = (-1)^h;
                Cond = [Cond, z(:,N+1) + d(i,:) * v - z(:,1:N) * lambda_feas(:,i) == 0];
            end
        else
            d = sdpvar(1);
            lambda_feas = sdpvar(N,1, 'full');
            Cond = [Cond, z(:,N+1)-z(:,1:N)*lambda_feas==0];
        end
        
        %
        for i= 1:N
            Cond = [ismember(x(:,i),X), ismember(z(:,i),Y), Cond];   
        end
        Cond = [ismember(x(:,N+1),X), Cond]; 
        %Control constraints
        for i= 1:N
            Cond = [ismember(u(:,i),U), Cond];
        end
        % System dynamics
        for i= 1:N
            Cond = [Cond,  x(:,i+1) - [A,B] * [x(:,i); u(:,i)] == 0, z(:,i) - C * x(:,i)==0];
        end
        Cond = [Cond, z(:,N+1) - C * x(:,N+1)==0];
        if strict ==1
            Cond = [Cond, sum(lambda_feas) == 1, d - exp(a) ==0 , lambda_feas(:) >= 0 ];
        else
            Cond = [Cond, sum(lambda_feas) == 1, lambda_feas - d >= 0, d >= 0];
        end
        % Lambda_vol = sdpvar(2*n, m, 'full');
        % Lambda_vol = sparse_stochastic(N, n_x, (n+1)/N);
        Cond_vol = [];
        % for i =1:n_x
        %     Cond_vol = [Cond_vol, D_1 * (V_X(:,i) - x(:,1:N) * Lambda_vol(:,i)) - r <= 0 ]; 
        % end

        Lambda_vol_y = sparse_stochastic(N, n_y, (m+1)/N);
        for i =1:n_y
            Cond_vol = [Cond_vol, D * (V_Y(:,i) - z(:,1:N) * Lambda_vol_y(:,i)) - r_2 <= 0 ]; 
        end
        % beta_vol = sdpvar(n, 1, 'full');
        % lambda = sdpvar(N,1,'full');
        % r = sdpvar(1);
        % Cond_vol = [Lambda_vol * X.H(:,1:end-1) == [eye(n);-eye(n)]/N];
        % for i=1:N
        %     Cond_vol = [Lambda_vol * X.H(:,end) <= r * lambda(i,:) * ones(2*n,1) + [eye(n);-eye(n)] * beta_vol , Cond_vol];
        % end
        % z = 0;
        % for j = 1:N
        %     z = z + lambda(i,:) * x(:,i) - beta_vol;
        % end
        % Cond_vol = [Cond_vol, sum(lambda) == 1, z == 0, lambda >= 0, Lambda_vol(:) >= 0, r >= 0];
        
        if strict ==1
            Obj =   norm(r_2, 2) - sum(a);
            prob = optimize([Cond_vol, Cond],Obj,ops);
        else
            Obj =   norm(r_2, 2) - d;
            prob = optimize([Cond_vol, Cond],Obj,ops);
        end
        
        if prob.problem ==0
            test = 1;
            x_1 = value(x);
            z_1 = value(z);
            fprintf("Min Conv_const: %.3f, Volume_output_space: %.3f \n", [min(value(d)), value(norm(r_2,2))])
            S_1 = update_vertices(z_1(:,1:N),x_1(:,1:N),V_null, A,B, C,U,X);
            H_poly = [H_poly S_1 x_1(:,1:N)]; 
            Z_poly = [Z_poly z_1(:,1:N)];
            %S = convhulln(H_poly');
            %fprintf("Number of rows: %.0f, Number of Column: %.0f , Number vertice: %.0f , Dimension: %.0f \n", [size(S), N,n]);
            H_ = Polyhedron('V', H_poly');
            H_.minVRep();
            H_poly = H_.V';
            if ~isequal(H_poly,[])
                c = c+1;
            end
        end
        if c == n_s || (~isequal(H_poly,[]) && N == N_max)
            %disp(H_poly);
            RCI = Polyhedron('V',transpose(H_poly));
            if RCI.isFullDim()
                disp("Initial Polyhedron is full dimensional");
            else
                disp("Initial Polyhedron is not full dimensional");
            end
            if method == 1 
                try
                    RCI.minHRep();
                    RCI = backward_compute(RCI,X, U, A, B, n_iter,n_x, N_max);
                catch error
                    disp(error.stack);
                    RCI = backward_compute(RCI,X, U, A, B, n_iter,n_x, N_max);
                end
                
            elseif method==2
                RCI = backward_compute_2(RCI,X, U, A, B, n_iter,1,L_check);
            end
            test = 1 ;
            break
        end
    end
    if test==0  
        warning("increase N_max value.");
        RCI =[];
    end
end
function Vec = update_vertices(Z_poly,H_poly,V_null, A,B, C,U,X)
    m = size(H_poly,2); 
    n = size(V_null,2); 
    f_star = 2*n*(m);
    Vec = zeros(size(H_poly,1),f_star);
    for i = 1:m
        h = {}; 
        u = {};
        v = {};
        alpha = {};
        z = {};
        Cond = {};
        lambda = {};
        c = {};
        h_star = zeros(1,m);
        z_vec = num2cell(repmat(H_poly(:,i), 1, 2*n),1);
        A_X = mat2cell(repmat(X.A, 1, 2*n),size(X.A,1),size(X.A,2)*ones(1,2*n));
        A_U = mat2cell(repmat(U.A, 1, 2*n),size(U.A,1),size(U.A,2)*ones(1,2*n));
        b_X = num2cell(repmat(X.b, 1, 2*n),1);
        b_U = num2cell(repmat(U.b, 1, 2*n),1);
        z_null = num2cell([V_null, -V_null],1);
        parfor j =1:2*n
            v{j} = z_null{j};
            u{j} =  sdpvar(size(B,2), 1);
            alpha{j} = sdpvar(1);
            z{j} = sdpvar(size(A,2), 1);
            x{j} = sdpvar(size(A,2), 1);
            b{j} = sdpvar(size(C,1), 1);
            c{j} = sdpvar(size(C,1), 1);
            lambda{j} = sdpvar(size(Z_poly,2), 1);
            lambda_1{j} = sdpvar(size(Z_poly,2), 1);
            Cond{j} = [A_X{j}*z{j} <= b_X{j}, A_X{j}*x{j} <= b_X{j}, A_U{j}*u{j} <= b_U{j}];
            Cond{j} = [Cond{j}, z{j} -  z_vec{j} - alpha{j} * v{j} == 0, x{j} -  [A,B] * [z{j}; u{j}] == 0];
            Cond{j} = [Cond{j}, C * z{j} - Z_poly * lambda{j} == 0,  C * x{j} - Z_poly * lambda_1{j} == 0];
            
            Cond{j} = [Cond{j}, alpha{j} >= 0];
            Obj{j} =  exp(-alpha{j});
            %ops{j} = sdpsettings('solver','ipopt','mu_strategy', 'adaptive', 'verbose',0);
            ops{j} = sdpsettings('solver','fmincon', 'verbose', 0);
            prob = optimize(Cond{j},Obj{j},ops{j});
            if prob.problem == 0
                z_vec{j} = value(z{j});
                alpha{j} = value(alpha{j});
            else 
                disp(prob.problem);
            end
            % [z_vec{j}, alpha{j}] = make_LP(v{j},z_vec{j}, A,B,C, Z_poly, A_X{j}, b_X{j}, A_U{j}, b_U{j})
        end
        %disp(max(alpha))
        h_star = 2*n*(i-1);
        Vec(:, h_star+1:h_star+2*n) = cell2mat(z_vec);
    end
end
function S = sparse_stochastic(n,m,density)
    S = round(sprand(n, m, density),2); % Sparse random matrix
    % Normalize columns to make it column stochastic
    colSums = sum(S, 1); % Sum of each column
    S = S * spdiags(1 ./ colSums', 0, m, m);
end
function x_opt = solve_LP(Cost, A, b)
    ops = optimoptions(@linprog,'display','final');
    [x,fval] = linprog(Cost,A,b,[],[],[],[], ops);
    x_opt = fval;
 end
 function l = distinct_randi(n, m)
    l = zeros([1,m], "int8");
    for i = 1:m
        h = randi(n,1);
        if ismember(l,h)
            while ismember(l,h)
                h = randi(n,1);
            end
        end
        l(i)=h;
    end
 end
 function [z_vec, alpha] = make_LP(v,z_vec, A,B,C, Z_poly, A_X, b_X, A_U, b_U)
    % Dimensions
    nz = size(A,2);      % z dimension
    nx = size(A,2);      % x dimension
    nu = size(B,2);
    nb = size(C,1);
    nc = size(C,1);
    nl = size(Z_poly,2);
    
    % Decision: y = [z; x; u; b; c; lambda; lambda1; alpha]
    idx_z      = 1:nz;
    idx_x      = nz + (1:nx);
    idx_u      = nz+nx + (1:nu);
    idx_b      = nz+nx+nu + (1:nb);
    idx_c      = nz+nx+nu+nb + (1:nc);
    idx_lam    = nz+nx+nu+nb+nc + (1:nl);
    idx_lam1   = nz+nx+nu+nb+nc+nl + (1:nl);
    idx_alpha  = nz+nx+nu+nb+nc+nl+nl + 1;
    
    Nvar = idx_alpha;
    
    %% Objective: min -alpha
    f = zeros(Nvar,1);
    f(idx_alpha) = -1;
    
    %% Equality constraints Aeq*y = beq
    Aeq = [];
    beq = [];
    
    % (1) b - Z_poly*lambda = 0
    tmp = zeros(nb,Nvar);
    tmp(:,idx_b)   = eye(nb);
    tmp(:,idx_lam) = -Z_poly;
    Aeq = [Aeq; tmp];
    beq = [beq; zeros(nb,1)];
    
    % (2) c - Z_poly*lambda1 = 0
    tmp = zeros(nc,Nvar);
    tmp(:,idx_c)    = eye(nc);
    tmp(:,idx_lam1) = -Z_poly;
    Aeq = [Aeq; tmp];
    beq = [beq; zeros(nc,1)];
    
    % (3) b - C*z = 0
    tmp = zeros(nb,Nvar);
    tmp(:,idx_b) = eye(nb);
    tmp(:,idx_z) = -C;
    Aeq = [Aeq; tmp];
    beq = [beq; zeros(nb,1)];
    
    % (4) c - C*x = 0
    tmp = zeros(nc,Nvar);
    tmp(:,idx_c) = eye(nc);
    tmp(:,idx_x) = -C;
    Aeq = [Aeq; tmp];
    beq = [beq; zeros(nc,1)];
    
    % (5) z - z_vec - alpha*v = 0
    tmp = zeros(nz,Nvar);
    tmp(:,idx_z) = eye(nz);
    tmp(:,idx_alpha) = -v;   % v is nz x 1
    Aeq = [Aeq; tmp];
    beq = [beq; z_vec];
    
    % (6) x - A*z - B*u = 0
    tmp = zeros(nx,Nvar);
    tmp(:,idx_x) = eye(nx);
    tmp(:,idx_z) = -A;
    tmp(:,idx_u) = -B;
    Aeq = [Aeq; tmp];
    beq = [beq; zeros(nx,1)];
    
    %% Inequality constraints Aineq*y <= bineq
    Aineq = [];
    bineq = [];
    
    % z in X: A_X*z <= b_X
    tmp = zeros(size(A_X,1),Nvar);
    tmp(:,idx_z) = A_X;
    Aineq = [Aineq; tmp];
    bineq = [bineq; b_X];
    
    % x in X: A_X*x <= b_X
    tmp = zeros(size(A_X,1),Nvar);
    tmp(:,idx_x) = A_X;
    Aineq = [Aineq; tmp];
    bineq = [bineq; b_X];
    
    % u in U: A_U*u <= b_U
    tmp = zeros(size(A_U,1),Nvar);
    tmp(:,idx_u) = A_U;
    Aineq = [Aineq; tmp];
    bineq = [bineq; b_U];
    
    % alpha >= 0  -->  -alpha <= 0
    tmp = zeros(1,Nvar);
    tmp(1,idx_alpha) = -1;
    Aineq = [Aineq; tmp];
    bineq = [bineq; 0];
    
    %% Solve LP
    opts = optimoptions('linprog','Display','none');
    [yopt, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, [], [], opts);
    
    if exitflag == 1
        z_vec = yopt(idx_z);
        alpha = yopt(idx_alpha);
    else
        z_vec = z_null;
        alpha = 0;
        disp(exitflag);
    end
end