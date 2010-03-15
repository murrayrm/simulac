function val = fate(data, MOI, PREmin)

% Hardcoded values for indexing
sl_time = 1;
sl_CI2 = 8+1;
sl_Cro2 = 9+1;
PREcols = [30, 32, 34, 36, 38, 40, 42, 44, 46, 48];

% Initialization
activeTwoMin = 0;			% PRE active in last two minutes?
lastActive = data(1, sl_time);		% Last time we saw an active PRE

% Figure out what the time step per iteration is
timeStep = data(2, sl_time) - data(1, sl_time);

% First condition: check if PRE is sufficiently activated
% Look through the simulation record step by step
for i = 1:size(data, 1)
  curTime = data(i, sl_time);

  % Make sure we have been running long enough
  if (i * timeStep < 120)
    % Not enough time elapsed to average
    continue
  end
        
  %
  % Compute the average number of RNAPs @ PRE over the last 2 min
  %

  % Figure out the beginning time for a two minute window
  startIndex = round(i - 120/timeStep) + 1;
  PREactive = mean(sum(data(startIndex:i, PREcols) == 3, 1));

  % See if we should start a timer for the 4 minute clock
  if ((PREactive > PREmin) && (activeTwoMin == 0))
    startFourMin = curTime;
    activeTwoMin = 1;
  end

  % Keep track of whether we have had active complex in last 2 min
  if (~(PREactive > 1))
    activeTwoMin = 0;
  end

  % See whether there has been an active complex for four minutes
  if (activeTwoMin && (curTime - startFourMin > 240))
    % Successful activation of PRE; now check [CI2) > [Cro2)
    if (data(size(data,1), sl_CI2) > data(size(data,1), sl_Cro2))
      val = 1; return;			% Lysogeny
    else
      val = 0; return; 			% Lysis
    end
  end
end

% If we get here, we failed to lysogenize
val = 0; return;
