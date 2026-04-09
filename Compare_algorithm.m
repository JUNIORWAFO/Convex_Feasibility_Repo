rng(42);
% Original system:
addpath("./Utils")
A = [0.5 1; 0 -0.5]; B = [1; 0];
% State constraints:
Gx =[   0.9147   -0.5402;
        0.2005    0.6213;
       -0.8193    0.9769;
       -0.4895   -0.8200;
        0.7171   -0.3581;
        0.8221    0.0228;
        0.3993   -0.8788];
Fx = [  0.5566;
        0.8300;
        0.7890;
        0.3178;
        0.4522;
        0.7522;
        0.1099];
D = Polyhedron('H',[Gx Fx]);
% Input constraints:
umin = -1;
umax = 1;
Gu = [1; -1]; Fu = [umax; -umin];
U = Polyhedron('H',[Gu Fu]);
% try 
%     tic
%     cis = computeRCIS(A,B,Gx,Fx,Gu,Fu,[],[],[],0,6);
%     toc
% catch
%     cis = Polyhedron('V', [0;0]);
% end
cis = Polyhedron('V', [0;0]);
%% Compute MCIS using MPT3: 
system = LTISystem('A',A,'B',B);
tic
mcis = system.invariantSet('X',D,'U',U,'maxIterations',300);
toc
%% Compute_using_feasibility
tic
CI_feas = compute_CI_feas(A,B, D,U,30,7,1,23, 2, 2,0, 1, 'ipopt' );
toc
% C = [0,1];
% Y = C*D;
% Y.computeHRep();
% Y.computeVRep();
% tic
%    CI_feas_obs = compute_CI_feas_Obs(A,B,C, D,U,Y,30,3,7,23,2,2, 0, 1,'fmincon');
% toc
% disp(CI_feas_obs<=mcis)
%% Plotting:
figure; plot(D, 'color', 'blue', cis(1), 'color', 'white', mcis, 'color', 'lightgray',CI_feas, 'color', 'red')
