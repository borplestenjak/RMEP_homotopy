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
diary Example_Arma2.txt

% set combinations of 
%   - m
%   - n
%   - number of runs of homotopy
%   - number of runs of multipareig_rect
%   - number of runs of macaulaylab
%   - estimated time (just for info for countdown)

Params_k_n = [
   1 4 1 1 0 10
   % 1 12 3 3 0 40
   % 1 4 5 5 5 20
   % 1 6 5 5 1 200
   % 1 8 5 5 0 25
   % 1 10 5 5 0 40
   % 1 14 3 3 0 65
   % 1 16 3 3 0 110
   % 1 18 3 3 0 200
   % 1 20 3 3 0 260
   % 1 22 1 1 0 130
   % 1 24 1 1 0 180
   % 1 26 1 1 0 250
   % 1 28 1 1 0 350
   % 1 30 1 1 0 460
   % 1 32 1 1 0 700
   % 1 34 1 1 0 1000
   % 2 6 3 3 0 50
   % 2 7 3 3 0 100
   % 2 8 3 3 0 200
   % 2 9 1 1 0 120
   % 2 10 1 1 0 230
   % 2 11 1 1 0 400
   % 2 12 1 1 0 800
   2 14 1 1 0 1600
   2 15 1 1 0 2500
   % 1 40 1 1 0 3000
   % 1 45 1 1 0 3000
   % 1 50 1 1 0 3000
    ];

Results = Params_k_n;
Results(end,10) = 0;

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
    [Res1, Res2, Res3] = compare_Arma_solvers(m,n,run1,run2,run3)
    if run1>0
        Results(j,8) = Res1.full_time_avg;
    end
    if run2>0
        Results(j,9) = Res2.full_time_avg;
    end
    if run3>0
        Results(j,10) = Res3.full_time_avg;
    end
    % we save results after every run if something goes wrong
    save Example_Arma2.mat Results
end

diary off