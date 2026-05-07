% runs big comparision of methods for linear RMEPS

% we save following results
% for all:
% - n, k, number of solutions
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
clc
diary Example_CS_Lin_RMEP.txt

% set combinations of 
%   - k
%   - n
%   - number of runs of homotopy
%   - number of runs of multipareig_rect
%   - number of runs of macaulaylab
%   - estimated time (just for info for countdown)

Params_k_n = [
    2  5 5 5 5 1;
    2 10 5 5 5 2;
    2 15 5 5 5 5;
    2 20 5 5 5 10;
    2 25 5 5 5 25;
    2 30 5 5 5 40;
    2 35 5 5 5 90;
    2 40 3 3 3 140;
    2 45 3 3 3 205;
    2 50 3 3 3 185;
    3 5 5 5 5 3; 
    3 10 5 5 5 10; 
    3 15 5 5 5 15; 
    3 20 3 3 3 375; 
    3 25 3 3 1 400; 
    3 30 3 3 0 390; 
    3 35 1 1 0 375; 
    ];

Results = Params_k_n;
Results(end,15) = 0;

for j = 1:size(Params_k_n,1)
    k = Params_k_n(j,1);
    n = Params_k_n(j,2);
    est_time = sum(Params_k_n(j:end,6));
    fprintf('Running example %d/%d with k=%d and n=%d, estimate time %d, estimate remaining time %d s\n',j,size(Params_k_n,1),k,n,Params_k_n(j,6),est_time)
    fprintf('========================================================================================\n')
    nsol = nchoosek(n+k-1,k);
    run1 = Params_k_n(j,3);
    run2 = Params_k_n(j,4);
    run3 = Params_k_n(j,5);
    Results(j,6) = nsol;
    [Res1, Res2, Res3] = compare_rmep_solvers(n,k,run1,run2,run3)
    if run1>0
        Results(j,7) = Res1.repeat_avg/nsol;
        Results(j,8) = Res1.avgsteps_avg;
        Results(j,9) = Res1.init_time_avg;
        Results(j,10) = Res1.full_time_avg;
        Results(j,11) = Res1.maxres_avg;
        Results(j,12) = Res1.meanres_avg;
    end
    if run2>0
        Results(j,13) = Res2.full_time_avg;
        Results(j,14) = Res2.maxres_avg;
        Results(j,15) = Res2.meanres_avg;
    end
    if run3>0
        Results(j,16) = Res3.full_time_avg;
        Results(j,17) = Res3.maxres_avg;
        Results(j,18) = Res3.meanres_avg;
    end
    % we save results after every run if something goes wrong
    save Example_CS_Lin_RMEP.mat Results
end

diary off