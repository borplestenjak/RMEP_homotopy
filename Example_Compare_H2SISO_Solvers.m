% runs big comparision of methods for polynomial RMEPS

% we save following results
% for all:
% - deg, n, k, number of solutions
% 1) homotopy
%       - number of repeated paths (in percentage)
%       - average number of steps
%       - time fot the construction of the initial problem
%       - overall time
%       - error
% 2) MultiparEig
%       - overall time
%       - error
% 3) Macaulay
%       - overall time
%       - error

% we use following combinations
% k, n, runs1, runs2, runs3
%
% Run parpool before this command

% Bor Plestenjak 2026

clear all
% clc
diary Example_H2SISO.txt

% set combinations of 
%   - m
%   - n
%   - number of runs of homotopy
%   - number of runs of multipareig_rect
%   - number of runs of macaulaylab
%   - estimated time (just for info for countdown)

Params_k_n = [
   2 5 1 0 0 10
   2 10 5 5 5 20
   2 30 3 3 0 200
   3 5 5 5 5 20
   3 10 3 3 3 200
   2 50 1 1 0 600
   3 15 1 1 1 2000
   4 10 1 0 1 1000
   4 15 1 0 0 1000
   2 30 0 0 1 1000
   2 70 1 1 0 3000
    ];

Results = Params_k_n;
Results(end,18) = 0;

for j = 1:size(Params_k_n,1)
    m = Params_k_n(j,1);
    n = Params_k_n(j,2);
    est_time = sum(Params_k_n(j:end,6));
    fprintf('Running example %d/%d with m=%d and n=%d, estimate time %d, estimate remaining time %d s\n',j,size(Params_k_n,1),m,n,Params_k_n(j,6),est_time)
    fprintf('========================================================================================\n')
    run1 = Params_k_n(j,3);
    run2 = Params_k_n(j,4);
    run3 = Params_k_n(j,5);
    % Results(j,7) = nsol;
    [Res1, Res2, Res3] = compare_H2SISO_solvers(m,n,run1,run2,run3)
    if run1>0
        Results(j,8) = Res1.full_time_avg;
        Results(j,9) = Res1.maxres_avg;
    end
    if run2>0
        Results(j,10) = Res2.full_time_avg;
        Results(j,11) = Res2.maxres_avg;
    end
    if run3>0
        Results(j,12) = Res3.full_time_avg;
        Results(j,13) = Res3.maxres_avg;
    end
    % we save results after every run if something goes wrong
    save Example_H2SISO.mat Results
end

diary off