clear all; clc;

ist = 1;
iat = 1.25;

logname = "BatchLog";
modelname = "MM1K";

W = 1000;       % Warm-up period
T = [10000];    % Batch size
C = 9;         % Number of bins
N = 50;         % Number of batches
M = 100;        % Samplepoints per batch for determining distribution in batch

global interServiceTime;
global interArrivalTime;
global stopTime;
global mySeed;
global capacity;

interServiceTime = Simulink.Parameter;
interArrivalTime = Simulink.Parameter;
stopTime = Simulink.Parameter;
mySeed = Simulink.Parameter;
capacity = Simulink.Parameter;

interServiceTime.Value = ist;
interArrivalTime.Value = iat;
capacity.Value = C;
% Determine expected pdf and cdf (the latter is needed for the KS test)
[expectedpdf, buffervalues] = pdfMM1K(ist/iat,C);
expectedcdf = pdf2cdf(expectedpdf);

starti = 1;
% In case the script is interrupted halfway execution, it can be continued by restoring the previous state at this point.
% Simply uncomment the next two lines in order to load all relevant variables and start the loop where it was interrupted.

% load(logname)
% starti = i+1;

% Start with an empty simulated pdf
simulatedcount = [];

% Running the simulations...
for i = starti : size(T,2)
    % Setting the simulation time
    stopTime = W + N * T(i);
    % Set a unique seed for the upcoming run of the simulator
    mySeed.Value = 0 + 2*i;
    % Run the simulator
    SimOut = sim(modelname,'ReturnWorkspaceOutputs', 'on');
    thisrun = SimOut.logsout{1}.Values;
    % Resample the run as a whole
    sampletimes = [W : T(i) / (M - 1) : W + N * T(i)];
    thisrun = resample(thisrun,sampletimes);
    % Split into batches and determine counts for each batch
    batchpoints = [W : T(i) : W + N * T(i)];
    batchpdf = [];
    for j = 1:size(batchpoints,2)-1
        batchruns(j) = getsampleusingtime(thisrun,batchpoints(j),batchpoints(j+1));
        batchpdf = [batchpdf ; sum([batchruns(j).Data == [0 : C - 1], batchruns(j).Data >= C])];
    end
    % Determine the mean of each count over all batches as the estimate of the pdf
    simulatedcount  = [simulatedcount ; mean(batchpdf)];
    % Plot the pdf by means of progress indication
    figure
    bar(buffervalues, simulatedcount(i,:),"blue");
    hold on;
    plot(buffervalues, expectedpdf * M,"red");
    title("M/M/1/" + (capacity.Value+1) + ...
        " Batch (ist =" + interServiceTime.Value + ...
        ", iat =" + interArrivalTime.Value + ...
        ", W = " + W + ", T = " + T(i) + ", N =" + N +", M = " + M + " )");
    xlabel("Number of tokens in buffer (expected = red, simulated = blue)");
    ylabel("Average number of occurences in batch");
    drawnow;
    save(logname,"i","j","interArrivalTime","interServiceTime","T","stopTime", ...
        "mySeed","simulatedcount","expectedpdf","C","N", "M", "W");
end