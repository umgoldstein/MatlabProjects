function GoldsteinProject1


% This function analysis the data from SCHATZ-study
%
%       Input: None
%
%       Ouptut: Answers to Project 1
%
%
% Project 1:
%
% Based on the entries in the runlog variable, the project task is to derive
% by means of Matlab programming the following list of experimental and behavioral indices,
% each in the form of an average over all participants, a standard deviation over participants, and the individual subject indices.
%
% - The average number of solvable tasks, i.e., tasks for which the number of possible steps allowed for the discovery of both treasures
% - The average task and run durations in
%        (a) number of decisions and (b) seconds
% - The average proportion of solved tasks across all solvable tasks
% - The average proportion of solved tasks across all solvable tasks that were solved with
%        (a) the optimal number of steps, (b) the optimal number of steps plus 1 step, (c) the optimal number of steps plus 2 steps
% - The median reaction time across all decisions
% - The median reaction time on attempts that followed an attempt on which a single treasure was discovered
% - Three further indices of your choice and design:
%   - the mean reaction time across all decision
%   - the average number of tasks solved on the frist attempt
%   - the average number of non-responses
%
% Copyright (C) Ursula Goldstein
% -------------------------------------------------------------------------
% clear up work space

clc
clear all

% define directory and data-----------------------------------------------

 src_dir = ('C:\Users\Ursula\Documents\MATLAB\Poject_1data\Matlab\Data')    ; % defined directory

 subj_dir = {'SCHATZ003',...
             'SCHATZ004',...
             'SCHATZ005',...
             'SCHATZ006',...
             'SCHATZ007',...
             'SCHATZ008',...
             'SCHATZ009',...
             'SCHATZ010',...
             'SCHATZ011',...
             'SCHATZ012',...
             'SCHATZ013',...
             'SCHATZ014',...
             'SCHATZ015',...
             'SCHATZ016',...
             'SCHATZ017',...
             'SCHATZ018',...
             'SCHATZ019',...
             'SCHATZ020',...
             'SCHATZ021',...
             'SCHATZ022'}                                                   ; % participants data file

 ses_num = { [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3 4],...
             [1 2 3],...
             [1 2 3 4] }                                                    ; % number of runs for each participant

 % Define empty vectors for later use ------------------------------------

           PRT = []                                                         ; % Participants median reaction time
        PR_std = []                                                         ; % Participants std for median reaction time
          PRT1 = []                                                         ; % Participants median reaction time after one treasure discovered
          PRT1 = []                                                         ; % Participants std (after one treasure discovered)
       PSolved = NaN(20,1)                                                  ; % Participants solved tasks
     PSolvable = NaN(20,1)                                                  ; % How many tasks were solvable in each task
      POptimal = NaN(20,1)                                                  ; % How many tasks were solved with optimal steps
     POptimal1 = NaN(20,1)                                                  ; % How many tasks were solved with optimal steps + 1
     POptimal2 = NaN(20,1)                                                  ; % How many tasks were solved with optimal steps + 2
  PRunAverageS = []                                                         ; % Average time taken for runs in seconds
 PTaskAverageS = NaN(20,1)                                                  ; % Average time taken for tasks in seconds
  PRunAverageD = []                                                         ; % Avergae number of decisions for runs
 PTaskAverageD = NaN(20,1)                                                  ; % Average number of decisions for tasks
       PRTmean = []                                                         ; % Participans mean reaction time
   PRTmean_std = []                                                         ; % Particiapnts std for mean reaction time
 PDontResponds = NaN(20,1)                                                  ; % How often did a participant not respond
   PSolved1st = NaN(20,1)                                                   ; % How many tasks were solved on the first attempt

% ------------------------------------------------------------------------
% ------------------------------ANALYSIS----------------------------------
% ------------------------------------------------------------------------

% funtcions for behavioral analysis
% ------------------------------------------------------------------------

for k = 1:numel(subj_dir)

    [tasks_solvable] = Solvable(src_dir, subj_dir{k},ses_num{k})            ; % function counts how many tasks were solvable
    PSolvable(k,1:length(tasks_solvable)) = tasks_solvable                  ; % collect data for solvable tasks for all participants

    [decisionrunaverage...
    ,decisiontaskaverage] = Decision(src_dir, subj_dir{k},ses_num{k})       ; % function counts the average number of decision in all tasks and runs
    PRunAverageD = [PRunAverageD;decisionrunaverage]                        ; % collect average run data for all participants
    PTaskAverageD(k,1:length(decisiontaskaverage)) = decisiontaskaverage    ; % collect average task data for all participants

    [taskaverage,runaverage] = Seconds(src_dir, subj_dir{k},ses_num{k})     ; % function counts the average time taken for runs and tasks
    PRunAverageS = [PRunAverageS;runaverage]                                ; % collect average run data for all participants
    PTaskAverageS(k,1:length(taskaverage)) = taskaverage                    ; % collect average task data for all participants

    [Solved12] = Solved(src_dir, subj_dir{k},ses_num{k})                    ; % function counts all tasks that were solved
 	PSolved(k,1:length(Solved12)) = Solved12                                ; % collect data for solved tasks for all participants

    [Optimal,Optimal1,Optimal2] = Optimal123(src_dir...
                                            ,subj_dir{k},ses_num{k})        ; % function counts how often solved with optimal steps or +1, +2
      POptimal(k,1:length(Optimal)) = Optimal                               ; % collect data for optimal steps (for all participants)
    POptimal1(k,1:length(Optimal1)) = Optimal1                              ; % collect data for optimal steps+1 "
    POptimal2(k,1:length(Optimal2)) = Optimal2                              ; % collect data for optimap steps+2 "

    [median_RT,RT_std] = medianRT(src_dir, subj_dir{k},ses_num{k})          ; % function for median reaction time for all participants
 	PRT(k,1) = median_RT                                                    ; % collect data for median reaction time for all participants
    PRT_std(k,1) = RT_std                                                   ; % collect data for median reaction time for all participants

	[medianRT1,RT1_std] = medianRTafter1(src_dir, subj_dir{k},ses_num{k})   ; % function for median reaction time when one treasure was previously discovered
    PRT1(k,1) = medianRT1                                                   ; % collect data for median reaction time for all participants
    PRT1_std(k,1) = RT1_std                                                 ; % collect data for median reaction time for all participants

	[mean_RT, meanRT_std] = meanRT(src_dir, subj_dir{k},ses_num{k})         ; % function for mean reaction time for all participants
 	PRTmean(k,1) = mean_RT                                                  ; % collect data for mean reaction time for all participants
    PRTmean_std(k,1) = meanRT_std                                           ; % collect data for mean reaction time for all participants

    [Solved1st] = Solved1stAttempt(src_dir, subj_dir{k},ses_num{k})         ; % function for how many tasks were solved on the 1st attempt
    PSolved1st(k,1:length(Solved1st)) = Solved1st                           ; % collect data for solved on 1st attempt for all participants

    [dontresponds] = DontRespond(src_dir, subj_dir{k},ses_num{k})           ; % function counts the non-responds
    PDontResponds(k,1:length(dontresponds)) = dontresponds                  ; % collect data in a vector from all participants

    [Tasks1] = CountTasks(src_dir, subj_dir{k},ses_num{k})                  ; % Function to count how many tasks there were
     Tasks(k,1) = Tasks1                                                    ; % collect data for CountTasks

end


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------
% Get mean and standard deviation over all participants ------------------

   [SolvedA] = SolvedAverage(subj_dir,PSolved,PSolvable)                    ; % function for average of solved tasks
  [OptimalA] = OptimalAverage(subj_dir,POptimal,PSolvable)                  ; % function for average of solved tasks with optimal steps
 [Optimal1A] = Optimal1Average(subj_dir,POptimal1,PSolvable)                ; % function for average of solved tasks with optimal steps +1
 [Optimal2A] = Optimal2Average(subj_dir,POptimal2,PSolvable)                ; % function for average of solved tasks with optimal steps +2

          PSolvable = (PSolvable./Tasks)*100                                ; % get Solvable to average over actual tasks
         PSolved1st = (PSolved1st./Tasks)*100                               ; % get Solved1st to average over actual tasks

      TaskSolvableA = mean(PSolvable)                                       ; % get mean for solvable tasks
    TaskSolvableStd = std(PSolvable)                                        ; % get standard deviation for solvable tasks

       DecisionRunA = mean(PRunAverageD)                                    ; % get mean for decision average for runs
     DecisionRunstd = std(PRunAverageD)                                     ; % get standard deviation for decision average for runs
      DecisionTaskA = mean(PTaskAverageD)                                   ; % get mean for decision average for tasks
    DecisionTaskstd = std(PTaskAverageD)                                    ; % get standard deviation for decision average for tasks

        SecondsRunA = mean(PRunAverageS)                                    ; % get mean for time in seconds for runs
      SecondsRunstd = std(PRunAverageS)                                     ; % get standard deviation for duration (s) for runs
       SecondsTaskA = mean(PTaskAverageS)                                   ; % get mean for time in seconds for tasks
     SecondsTaskstd = std(PTaskAverageS)                                    ; % get standard deviation for duration (s) for tasks

           SolvedAA = mean(SolvedA)                                         ; % get mean of averaged solved tasks
           SolvedAS = std(SolvedA)                                          ; % get standard deviation of averaged solved tasks

          OptimalAA = mean(OptimalA)                                        ; % get mean of averaged solved tasks with optimal steps
          OptimalAS = std(OptimalA)                                         ; % get standard deviation of averaged solved tasks with optimal steps

         Optimal1AA = mean(Optimal1A)                                       ; % get mean of averaged solved tasks with optimal steps +1
         Optimal1AS = std(Optimal1A)                                        ; % get standard deviation of averaged solved tasks with optimal steps +1

         Optimal2AA = mean(Optimal2A)                                       ; % get mean of averaged solved tasks with optimal steps +2
         Optimal2AS = std(Optimal2A)                                        ; % get standard deviation of averaged solved tasks with optimal steps +2

          MedianRTA = mean(PRT)                                             ; % get mean for median reaction time
        MedianRTStd = std(PRT)                                              ; % get standard deviation for median reaction time

         MedianRT1A = mean(PRT1)                                            ; % get mean for median reaction time after one treasure was discovered
       MedianRT1Std = std(PRT1)                                             ; % get standard deviation for median reaction time after one treasure was discovered

            MeanRTA = mean(PRTmean)                                         ; % get mean for mean reaction time
          MeanRTStd = std(PRTmean)                                          ; % get standard deviation for mean reaction time

         Solved1stA = mean(PSolved1st)                                      ; % get mean for solved on 1st attempt
       Solved1stStd = std(PSolved1st)                                       ; % get standard deviation for solved on 1st attempt

    DontRespondMean = mean(PDontResponds)                                   ; % get mean of non-responds
    DontRespondStd  = std(PDontResponds)                                    ; % get standard deviation of non-responds

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% clear up workspace from unnessecary variables --------------------------
% ------------------------------------------------------------------------

clear bsc decisionrunaverage decisiontaskaverage dontresponds k mean_RT median_RT...
    medianRT1 Optimal Optimal1 Optimal2 runaverage Solved12 Solved1st taskaverage tasks_solvable...
    ses_num src_dir subj_dir Tasks1 RT1_std RT_std meanRT_std

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ---------------------------PRINT ANSWERS--------------------------------
% ------------------------------------------------------------------------

printmat(PSolvable, 'The average number of solvable tasks,',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', '%')
fprintf('The average number of solvable tasks overall participants is %.2d with a standard deviation of %.2f\n\n\n',...
    TaskSolvableA,TaskSolvableStd)

printmat(PRunAverageD, 'The average run duration in number of decisions ',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'DecisionRun')
fprintf('The average run duration in number of decisions overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    DecisionRunA,DecisionRunstd)

printmat(PTaskAverageD, 'The average task duration in number of decisions',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'DecisionTask')
fprintf('The average task duration in number of decisions overall participant is %.2f with a standard deviation of %.2f\n\n\n',...
    DecisionTaskA,DecisionTaskstd)

printmat(PRunAverageS, 'The average run duration in seconds',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 's')
fprintf('The average run duration in seconds overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    SecondsRunA,SecondsRunstd)

printmat(PTaskAverageS, 'The average task duration in seconds',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 's')
fprintf('The average task duration in seconds overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    SecondsTaskA,SecondsTaskstd)

printmat(SolvedA, 'The average proportion of solved tasks across all solvable tasks',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', '%')
fprintf('The average proportion of solved tasks across all solvable tasks that were solved overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    SolvedAA,SolvedAS)

printmat(OptimalA, 'The average proportion of solved tasks across all solvable tasks that were solved with the optimal number of steps ',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'Optimal')
fprintf('The average proportion of solved tasks across all solvable tasks that were solved with the optimal number of steps overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    OptimalAA,OptimalAS)

printmat(Optimal1A, 'The average proportion of solved tasks across all solvable tasks that were solved with the optimal number of steps plus 1 step ',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'Optimal+1')
fprintf('The average proportion of solved tasks across all solvable tasks that were solved with the optimal number of steps plus 1 step overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    Optimal1AA,Optimal1AS)

printmat(Optimal2A, 'The average proportion of solved tasks across all solvable tasks that were solved with the optimal number of steps plus 2 step ',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'Optimal+2')
fprintf('The average proportion of solved tasks across all solvable tasks that were solved with the optimal number of steps plus 2 step overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    Optimal2AA,Optimal2AS)

printmat(PRT, 'The median reaction time across all decisions',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'ms')
fprintf('The median reaction time across all decisions overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    MedianRTA,MedianRTStd)

printmat(PRT1, 'The median reaction time on attempts that followed an attempt on which a single treasure was discovered',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'ms')
fprintf('The median reaction time on attempts that followed an attempt on which a single treasure was discovered overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    MedianRT1A,MedianRT1Std)

printmat(PRTmean, 'The mean reaction time across all decisions',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'ms')
fprintf('The mean reaction time across all decisions overall participants is %.2f with a standard deviation of %.2f\n\n\n',...
    MeanRTA,MeanRTStd)

printmat(PSolved1st, 'How often did a participant solve the task on the first attempt',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'Tasks')
fprintf('The average number participants did solve the task on the first attempt on average is %.2f with a standard deviation of %.2f\n\n\n',...
    Solved1stA,Solved1stStd)

printmat(PDontResponds, 'How often did participants not respond in all 4 runs',...
    '003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022', 'NotResponded')
fprintf('The average number participants did not respond on average is %.2f on average with a standard deviation of %.2f\n\n\n',...
    DontRespondMean,DontRespondStd)

% -----------------------------------------------------------------------
% -----------------------------Generate Graphs---------------------------
% -----------------------------------------------------------------------

figure('name','Run Duration on Average in Decisions','numbertitle','off')

hold on
bar(PRunAverageD)
ylabel('Decisions')
xlabel('Participants')
xlim([0 21])
hold off

figure('name','Task Duration on Average in Decisions','numbertitle','off')

hold on
bar(PTaskAverageD)
ylabel('Decisions')
xlabel('Participants')
xlim([0 21])
ylim([0 70])
hold off

figure('name','Run Duration on Average in Seconds','numbertitle','off')

hold on
bar(PRunAverageS)
ylabel('time (s)')
xlabel('Participants')
xlim([0 21])
hold off

figure('name','Task Duration on Average in Seconds','numbertitle','off')

hold on
bar(PTaskAverageS)
ylabel('time (s)')
xlabel('Participants')
xlim([0 21])
hold off

figure('name','Tasks Solvable on Average','numbertitle','off')

hold on
bar(PSolvable)
ylabel('%')
xlabel('Participants')
ylim([0 100])
xlim([0 21])
hold off

figure('name','Tasks Solved on Average','numbertitle','off')

hold on
bar(SolvedA)
ylabel('%')
xlabel('Participants')
ylim([0 100])
xlim([0 21])
hold off

figure('name','Solved with Optimal Steps','numbertitle','off')

hold on
bar(OptimalA)
ylabel('%')
xlabel('Participants')
ylim([0 100])
xlim([0 21])
hold off

figure('name','Solved with Optimal Steps+1','numbertitle','off')

hold on
bar(Optimal1A)
ylabel('%')
xlabel('Participants')
ylim([0 100])
xlim([0 21])
hold off

figure('name','Solved with Optimal Steps+2','numbertitle','off')

hold on
bar(Optimal2A)
ylabel('%')
xlabel('Participants')
ylim([0 100])
xlim([0 21])
hold off

figure('name','Median Reaction Time','numbertitle','off')

hold on
bar(PRT)
ylabel('RT (ms)')
xlabel('Participants')
ylim([0 1400])
xlim([0 21])
errorbar(PRT,PRT_std, 'kx')
hold off

figure('name','Median Reaction Time After 1 Treasure found','numbertitle','off')

hold on
bar(PRT1)
ylabel('RT (ms)')
xlabel('Participants')
ylim([0 1400])
xlim([0 21])
errorbar(PRT1,PRT1_std, 'kx')
hold off

figure('name','Mean Reation Time','numbertitle','off')

hold on
bar(PRTmean)
ylabel('RT (ms)')
xlabel('Participants')
xlim([0 21])
ylim([0 1400])
errorbar(PRTmean,PRTmean_std, 'kx')
hold off

figure('name','Solved on 1st Attempt','numbertitle','off')

hold on
bar(PSolved1st)
ylabel('%')
xlabel('Participants')
ylim([0 100])
xlim([0 21])
hold off

figure('name','Did not Respond','numbertitle','off')

hold on
bar(PDontResponds)
ylabel('Dont Responds')
xlabel('Participants')
ylim([0 15])
xlim([0 21])
hold off



end

% ------------------------------------------------------------------------
% ----------------------------------END-----------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
%
% Following functions are used in the behavioral analysis above
%
%

function [tasks_solvable] =  Solvable(src_dir, subj_dir,ses_num)

% this function calculates the average number of solvable tasks for a
% participant.
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT: tasks_solvable - a vector containing the number of
%                      solvable tasks

%--------------------------------------------------------------------------
%

% Define empty vectors

tasks_solvable1 = []                                                        ;

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

solvable1 = []                                                              ; % clear vector for each run

j = 1                                                                       ; % set j to zero for each run
for i = 1:length(runlog)
    trial = sum(~isnan(runlog{1,i}(:,8)))                                   ; % get number of trials
    a = runlog{1,i}(trial,10)                                               ; % important for data location
    if a == 0
       c = runlog{1,i}(trial,14)                                            ; % location for steps
       d = runlog{1,i}(trial,15)                                            ; % location for optimal steps
       e = runlog{1,i}(trial,12)                                            ; % location for when tasks ended
       if c>=d
           solvable1(j,i) = 1                                               ; % count if solvable
           if e == 1 || e == 3
               j = j+1                                                      ; % count if task ended
           end
       else
           solvable1(j,i) = 0                                               ; % count if solvable
           if e == 1 || e == 3
               j = j+1                                                      ; % count if task ended
           end
       end
    else
        c = runlog{1,i}(trial+1,14)                                         ; % location for steps
        d = runlog{1,i}(trial+1,15)                                         ; % location for optimal steps
        e = runlog{1,i}(trial+1,12)                                         ; % location for when tasks ended
        if c>=d
            solvable1(j,i) = 1                                              ; % count if solvable
            if e == 1 || e == 3
               j = j+1                                                      ; % count if task ended
            end
        else
            solvable1(j,i) = 0                                              ; % count if solvable
           if e == 1 || e == 3
               j = j+1                                                      ; % count if task ended
           end
        end
    end
end

tasks_solvable2 = 0                                                         ; % set tasks counter to zero

for k = 1:length(s_opt_all)

    f = sum(solvable1(k,:))                                                 ; % count how many were solvable
        if f >= 1
            tasks_solvable2 = tasks_solvable2+1                             ; % count how many were solvable
        end
end

tasks_solvable1(r) = tasks_solvable2                                        ; % collect data in vector

end
tasks_solvable = sum(tasks_solvable1)                                       ; % sum how many tasks were solvable
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [decisionrunaverage,...
          decisiontaskaverage]  =  Decision(src_dir, subj_dir,ses_num)

% this function calculates the average for the number of decisions in all
% tasks and runs for a participant
%
%           INPUT: src_dir, subj_dir, ses_sum
%
%          OUTPUT:  decisionrunaverage - a vector containing all data from
%                                        a participant for run average
%                  decisiontaskaverage - a vector containing all data from
%                                        a participant for task average
%
% Define empty vectors
% --------------------------------------------
       decisiontask = []                                                    ;
decisiontaskaverage = []                                                    ;
              tasks = 0                                                     ;


for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

 decisions = 0                                                              ; % set counter to zero before each participant

    for i = 1:length(runlog)
        trial = sum(~isnan(find((runlog{1,i}(:,12))>3)))                    ; % get number of trial
        decisions = sum(~isnan(runlog{1,i}(1:trial,16)))+decisions          ; % sum all cells that are numbers (made a decision)

    end

tasks = length(s_opt_all)+tasks                                             ; % count number of tasks
decisiontask(r) = decisions                                                 ; % fill in number decision into vector

end

decisiontaskaverage = (sum(decisiontask))/tasks                             ; % get average overall tasks
decisionrunaverage = (sum(decisiontask))/(length(ses_num))                  ; % get average overall participants

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [taskaverage,...
          runaverage] =  Seconds(src_dir, subj_dir,ses_num)

% this function calculates the duratin of a task and run in seconds
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT: taskaverage - a vector containing the average
%                                    duration in second for a task
%                       runaverage - a vector containing the average
%                                    duration in second for a run
%
%--------------------------------------------------------------------------
%

% Define empty vectors

taskaverage1 = []                                                           ;
  runaverage = 0                                                            ;
       tasks = 0                                                            ;

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

    task_duration = 0                                                       ; % set counter to zero for each run
    trial = sum(~isnan(runlog{1,(length(runlog))}(:,8)))                    ; % get number of trials
    a     = runlog{1,(length(runlog))}(trial,10)                            ; % important for location of data

        if a == 0
            task_duration = runlog{1,(length(runlog))}(trial,18)+task_duration; % get task duration
        else
            task_duration = runlog{1,(length(runlog))}(trial+1,18)+task_duration; % get task duration
        end

          tasks = length(s_opt_all)+tasks                                   ; % count number of tasks

taskaverage1(r) = task_duration                                             ; % collect data in vector
end
    taskaverage = ((sum(taskaverage1))/tasks)/1000                          ; % convert data from ms to s and average over tasks
     runaverage = ((sum(taskaverage1))/(length(ses_num)))/1000              ; % convert data from ms to s and average over runs
end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [Solved12] =  Solved(src_dir, subj_dir,ses_num)

% this function calculates how often a participant solved a task
%
%               INPUT: src_dir, subj_dir,ses_num
%
%              OUTPUT: Solved12 - a vector containing how often a
%                                 participant solved a task in a run
%
%--------------------------------------------------------------------------
%

% Define empty vectors

Solved1 = []                                                                ;

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

solved12 = []                                                               ; % clear vector for run

for i = 1:length(runlog)
    trial = sum(~isnan(runlog{1,i}(:,8)))                                   ; % get number of trials
    a = runlog{1,i}(trial,10)                                               ; % important for data location
    if a == 0
       e = runlog{1,i}(trial,12)                                            ; % see if tasks was solved
           if e == 1
               solved12(i) = 1                                              ; % tasks was solved
           else
               solved12(i) = 0                                              ; % task was not solved
           end
    else
        e = runlog{1,i}(trial+1,12)                                         ; % see if tasks was solved
            if e == 1
                solved12(i) = 1                                             ; % tasks was solved
            else
                solved12(i) = 0                                             ; % task was not solved
            end
    end
end


% taskcounter--------------------------------------------------------------
% counts how many tasks there were, and how many attempts

task_1 = []                                                                 ; % set task counter to zero

for i = 1:length(runlog)
    trial = sum(~isnan(runlog{1,i}(:,8)))                                   ; % get number of trials
    a = runlog{1,i}(trial,10)                                               ; % important for data location
    if a == 0
       e = runlog{1,i}(trial,12)                                            ; % when was the tasks ended
       if e==1 || e==3
           task_1(i) = 1                                                    ; % task was ended
       end
    else
        e = runlog{1,i}(trial+1,12)                                         ; % when was the tasks ended
        if e==1 || e==3
            task_1(i) = 1                                                   ; % task was ended
        end
    end
end
    A = task_1 + solved12                                                   ; %
    Solved123 = sum(sum(A == 2))                                            ; % see how many were solved
    Solved1(r) = Solved123                                                  ; % collect data

end
    Solved12 = sum(Solved1)                                                 ; % sum how many tasks were solved
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [Optimal,...
          Optimal1,...
          Optimal2] =  Optimal123(src_dir, subj_dir,ses_num)

% this function calculates how often a participant solved a task with
% optimal steps, optimal steps +1, or optimal steps+ 2
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT:  Optimal - vector containing how often a tasks was
%                                 solved with optimal steps
%                      Optimal1 - vector containing how often a tasks was
%                                 solved with optimal steps +1
%                      Optimal2 - vector containing how often a tasks was
%                                 solved with optimal steps +2
%
%--------------------------------------------------------------------------

% Define empty vectors

 Optima0l = []                                                              ;
Optimal11 = []                                                              ;
Optimal22 = []                                                              ;

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

   solved = 0                                                               ; % set all counter to zero for each run
Optimal_0 = 0                                                               ;
Optimal_1 = 0                                                               ;
Optimal_2 = 0                                                               ;


for i = 1:length(runlog)
    trial = sum(~isnan(runlog{1,i}(:,8)))                                   ; % get number of trials
    a = runlog{1,i}(trial,10)                                               ; % important for location of data
    steps_taken = trial-1                                                   ; % number of steps taken
    if a == 0
       d = runlog{1,i}(trial,15)                                            ; % location of optimal steps
       e = runlog{1,i}(trial,12)                                            ; % needed to see if solved
           if e == 1
               solved = solved+1;
               Optimal_0 = (steps_taken == d)+Optimal_0                     ; % if the same as d, then optimal steps were used
               Optimal_1 = (steps_taken-1 == d)+Optimal_1                   ; % if the same as d-1, then optmal steps+1 were used
               Optimal_2 = (steps_taken-2 == d)+Optimal_2                   ; % if the same as d-2, then optmal steps+2 were used
           end
    else
        d = runlog{1,i}(trial+1,15)                                         ; % location of optimal steps
        e = runlog{1,i}(trial+1,12)                                         ; % needed to see if solved
            if e == 1
                solved = solved+1;
                Optimal_0 = (steps_taken == d)+Optimal_0                    ; % if the same as d, then optimal steps were used
                Optimal_1 = (steps_taken-1 == d)+Optimal_1                  ; % if the same as d-1, then optmal steps+1 were used
                Optimal_2 = (steps_taken-2 == d)+Optimal_2                  ; % if the same as d-2, then optmal steps+2 were used
            end
    end
end

 Optima0l(r) = Optimal_0                                                    ; % collect data in a vector
Optimal11(r) = Optimal_1                                                    ; % collect data in a vector
Optimal22(r) = Optimal_2                                                    ; % collect data in a vector
end
     Optimal = sum(Optima0l)                                                ; % sum data in the vector
    Optimal1 = sum(Optimal11)                                               ; % sum data in the vector
    Optimal2 = sum(Optimal22)                                               ; % sum data in the vector
end

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function [median_RT,RT_std] =  medianRT(src_dir, subj_dir,ses_num)

% this function calculates the median reaction time for a participant
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT: median_RT - a vector containing the median
%                                  reaction time of a participant
%
%--------------------------------------------------------------------------

% Define empty vectors

median_RT2 = []                                                             ;
 median_RT = 0                                                              ;

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

median_RT1 = NaN(14,12)                                                     ; % define NaN-vecotr for each run

for i = 1:length(runlog)

     trial = sum(~isnan(find((runlog{1,i}(:,12))>3)))                       ; % get number of trials
    median_RT1(1:length(runlog{1,i}(1:trial,16)),i)= runlog{1,i}(1:trial,16); % collect reaction times

end
median_RT2 = [median_RT2,median_RT1]                                        ; % collect data in a vector
end

 median_RT = nanmedian(median_RT2(:))                                       ; % calculate median (without NaNs)
 RT_std = nanstd(median_RT2(:))                                             ; % calculate standard deviation (w/o NaNs)

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [medianRT1, RT1_std] =  medianRTafter1(src_dir, subj_dir,ses_num)

% this function calculates the median reaction time for a participant after
% one treasure was discovered
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT: medianRT1 - a vector containing the median
%                                  reaction time of a participant after one
%                                  treasure was discovered in previous
%                                  attempt
%
%--------------------------------------------------------------------------

% Define empty vectors

  medianRT1 = 0                                                             ;
 median_RT2 = []                                                            ;

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

 median_RT3 = NaN(14,12)                                                    ; % define NaN vector for each run

for i = 1:length(runlog)
    trial = sum(~isnan(runlog{1,i}(:,8)))                                   ; % get number of trials
    a = runlog{1,i}(trial,5)                                                ; % cell for treasure 1 found or not
    b = runlog{1,i}(trial,6)                                                ; % cell for treasure 2 found or not
    e = runlog{1,i}(trial,10)                                               ; % important for location of reaction time
    if a == 0 && b == 1 || a == 1 && b == 0
        if e == 0
            d = runlog{1,i}(trial,12)                                       ; % if there was another attempt
            if d == 2
                trial1 = sum(~isnan(find((runlog{1,i+1}(:,12))>3)))         ; % if so get median reaction time
                median_RT3(1:length(runlog{1,i+1}(1:trial1,16)),i)= runlog{1,i+1}(1:trial1,16);
            end
        else
            d = runlog{1,i}(trial+1,12)                                     ; % if there was another attempt
            if d == 2
                trial1 = sum(~isnan(find((runlog{1,i+1}(:,12))>3)))         ; % if so get median reaction time
                median_RT3(1:length(runlog{1,i+1}(1:trial1,16)),i)= runlog{1,i+1}(1:trial1,16);
            end
        end
    end
end
median_RT2 = [median_RT2,median_RT3]                                        ; % collect reaction time
end
medianRT1 = nanmedian(median_RT2(:))                                        ; % get median (without NaNs)
RT1_std = nanstd(median_RT2(:))                                             ; % get standard deviation (without NaNs)
end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [mean_RT,meanRT_std] =  meanRT(src_dir, subj_dir,ses_num)

% This function calculates the mean reaction time for a participant
%
%           INPUT: src_dir, subj_dir,ses_num
%
%          OUTPUT: mean_RT - vector containing mean reaction time for a
%                            participant
%
% ------------------------------------------------------------------------

% define empty vectors and counters
mean_RT2 = []                                                               ;
 mean_RT = 0                                                                ; % set counter to zero for each participant

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

mean_RT1 = NaN(14,12)                                                       ; % define NaN-vector for each run

for i = 1:length(runlog)

     trial = sum(~isnan(find((runlog{1,i}(:,12))>3)))                       ; % get number of trials
     mean_RT1(1:length(runlog{1,i}(1:trial,16)),i)= runlog{1,i}(1:trial,16) ; % collect reaction times

end
mean_RT2 = [mean_RT2,mean_RT1]                                              ; % collect data in a vector
end
 mean_RT = nanmean(mean_RT2(:))                                             ; % get mean (except NaNs) for vector
 meanRT_std = nanstd(mean_RT2(:))                                           ; % get std (w/o NaNs)
end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [Solved1st] =  Solved1stAttempt(src_dir, subj_dir,ses_num)

% this function calculates how often a participant solved a task on the
% first attempt
%
%               INPUT:  subj_dir, PSolved, PSolvable
%
%              OUTPUT: Solved1st - a vector containing how often a
%                                  participant solved a task on the 1st attempt
%
%--------------------------------------------------------------------------
%

% Define empty vectors

Solved1st3 = []                                                             ;
for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

solved1st2 = 0                                                              ; % set counter to zero for each run

for i = 1:length(runlog)
    trial = sum(~isnan(runlog{1,i}(:,8)))                                   ; % get number of trials
    a = runlog{1,i}(trial,10)                                               ; % important for location of data
    if a == 0
       e = runlog{1,i}(trial,12)                                            ; % important to see if task solved
       f = runlog{1,i}(1,9)                                                 ; % important to see if  1st attempt
           if e == 1 && f == 1
               solved1st2 = solved1st2+1                                    ; % count if solved on the 1st attempt
           end
    else
        e = runlog{1,i}(trial+1,12)                                         ; % important to see if task solved
        f = runlog{1,i}(1,9)                                                ; % important to see if  1st attempt
           if e == 1 && f == 1
               solved1st2 = solved1st2+1                                    ; % count if solved on the 1st attempt
           end
    end
end
Solved1st3(r) = solved1st2                                                  ; % collect data in a vector
end
Solved1st = sum(Solved1st3)                                                 ; % sum collected data
end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [dontresponds]  =  DontRespond(src_dir, subj_dir,ses_num)

% this function counts how many times a participant did not respond and did
% not make a decision
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT: dontresponds - a vector containing the non-responds
%                                     of a participant
%
%--------------------------------------------------------------------------



dontresponds = 0                                                            ; % set counter to zero for each participant
for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

    for i = 1:length(runlog)
        trial = sum(~isnan(find((runlog{1,i}(:,12))>3)))                    ; % get number of trials
        dontresponds = sum(isnan(runlog{1,i}(1:trial,16)))+dontresponds     ; % count non-responds by summing NaNs
    end


end
end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


function [ SolvedA ] = SolvedAverage(subj_dir,PSolved,PSolvable)

% this function calculates the average of solved tasks over solvable tasks
%
%               INPUT:  subj_dir, PSolved, PSolvable
%
%              OUTPUT: SolvedAverage - a vector containing the PSolved data
%                                      averaged over the PSolvable data
%
%--------------------------------------------------------------------------
%

% Define empty vectors

SolvedA = []                                                                ;

for i = 1:numel(subj_dir)
SolvedA = [SolvedA;PSolved(i)./PSolvable(i)]                                ; % average solved tasks over solvable
end

SolvedA = SolvedA*100                                                       ; % get to percentage

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [OptimalA] = OptimalAverage(subj_dir,POptimal,PSolvable)
% this function calculates how often a participant solved a task with
% optimal steps overall solvable tasks
%
%               INPUT:  subj_dir, POptimal, PSolvable
%
%              OUTPUT:  OptimalA - vector containing the Optimal data
%                                  averaged over solvable tasks
%
%--------------------------------------------------------------------------
%

OptimalA = []                                                               ; % Define empty vector

for i = 1:numel(subj_dir)
OptimalA = [OptimalA;POptimal(i)./PSolvable(i)]                             ; % devide through solvable tasks
end

OptimalA = OptimalA*100                                                     ; % get to percentage

end

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

function [Optimal1A] = Optimal1Average(subj_dir,POptimal1,PSolvable)

% this function calculates how often a participant solved a task with
%  optimal steps +1 averaged over solvable tasks
%
%               INPUT:  subj_dir, POptimal1, PSolvable
%
%              OUTPUT:  Optimal1A - vector containing the Optimal1 data
%                                  averaged over solvable tasks
%
%--------------------------------------------------------------------------


Optimal1A = []                                                              ; % Define empty vector

for i = 1:numel(subj_dir)
Optimal1A = [Optimal1A;POptimal1(i)./PSolvable(i)]                          ; % Average Optimal1 over solvable tasks
end

Optimal1A = Optimal1A*100                                                   ; % get to percentage

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


function [Optimal2A] = Optimal2Average(subj_dir,POptimal2,PSolvable)

% this function calculates how often a participant solved a task with
% optimal steps+ 2, averaged over solvable tasks
%
%               INPUT:  subj_dir, POptimal2, PSolvable
%
%              OUTPUT:  Optimal2A - vector containing the Optimal2 data
%                                   averaged over solvable tasks
%
%--------------------------------------------------------------------------
%

Optimal2A = []                                                              ; % Define empty vector

for i = 1:numel(subj_dir)
Optimal2A = [Optimal2A;POptimal2(i)./PSolvable(i)]                          ; % Average Optimal1 over solvable tasks
end

Optimal2A = Optimal2A*100                                                   ; % get to percentage

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [ Tasks1 ] = CountTasks(src_dir, subj_dir,ses_num)
% This function counts how many tasks there were
%
%               INPUT: src_dir, subj_dir, ses_num
%
%              OUTPUT: Tasks - vector containing number of tasks
%
%--------------------------------------------------------------------------


Tasks2 = [];

for r = 1:length(ses_num)

    h = fullfile(src_dir, subj_dir, [subj_dir, '_Run_' num2str(r), '.mat']) ; % load data

    load(h)

   Tasks2(r) = length(s_opt_all) ;

end

Tasks1 = sum(Tasks2);

end
