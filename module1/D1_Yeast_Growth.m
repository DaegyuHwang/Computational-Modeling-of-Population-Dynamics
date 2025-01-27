% Simulation of Yeast Growth model (Continuously growing species)

% Under optimal growth conditions, it can have a doubling time as fast as
% 90minutes => let tau2 = 1.5hours(90minutes)

% David Hwang, 09/16/24


%%%%%%%%%%%%%%%%%%%%
%% Unrestricted(exponential) growth model 
%%%%%%%%%%%%%%%%%%%%

% main parameters

N0 = 1; % Beginning with a single yeast cell
timev = 0:3:24*3; % time vector for every 3hours
tau2 =  1.5; % doubling time of yeast is 1.5hours
N = N0; % N is a number; N0 is a initial value 
r = 3/(1e+6); % radius of Yeast (m)

% Simulation
sim = zeros(length(timev),1); % vector to store the values of N during simulation
k=0; % counter

for t=timev
    k=k+1;

    % main equation
    N = N0*2^(t/tau2);

    % store value of N for plotting
    sim(k)=N;
end

sim(k) % the number of cells at the end of the simulation; k=25

vol = 4/3*pi*r^3; % volume of a Yeast(m^3)
L = vol*10^(3); % liter of a Yeast(L)
L*sim(k) % Liters occupied by yeast at the end of the simulation(L); k=25

% plot
figure(1); clf
plot(timev,sim,'-o')
xlabel('time (hours)');ylabel('population (cells)')
title('Unrestricted(exponential) growth model for Yeast')




%%%%%%%%%%%%%%%%%%%%
%% Unrestricted(exponential) growth model 
%%%%%%%%%%%%%%%%%%%%

% main parameters

N0 = 10000; % Beginning with a single yeast cell
timev = 0:1:24*5; % time vector for every one hours
tau2 =  1.5; % doubling time of yeast is 1.5hours
N = N0; % N is a number; N0 is a initial value 
r = 3/(1e+6); % radius of Yeast (m)
vol = 4/3*pi*r^3; % volume of a Yeast(m^3)
L = vol*10^(3); % liter of a Yeast(L)


% Simulation
k=0; % counter
max_t=0; % hours to fill 1L

for t=timev
    k=k+1;

    % main equation
    N = N0*2^(t/tau2);

    % find out how many hours it takes to fill 1L
    if L*N>=1 max_t=t;break; end    

end



%%%%%%%%%%%%%%%%%%%%
%%  considering competition
%%%%%%%%%%%%%%%%%%%%

max_d1 = 2*10^8; % maximal density of yeast cells per mL
max_d2 = max_d1*1000; % maximal density of yeast cells per L


r = 3/(1e+6); % radius of Yeast (m)
vol = 4/3*pi*r^3; % volume of a Yeast(m^3)
L = vol*10^(3); % liter of a Yeast(L)

L*max_d2
vol*max_d2


%%%%%%%%%%%%%%%%%%%%
%% Logistic equation
%%%%%%%%%%%%%%%%%%%%

% main parameters

N0 = 10000; % Beginning with a single yeast cell
timev = 0:1:24*3; % time vector for every one hours
tau2 =  1.5; % doubling time of yeast is 1.5hours
N = N0; % N is a number; N0 is a initial value 
r = 3/(1e+6); % radius of Yeast (m)
vol = 4/3*pi*r^3; % volume of a Yeast(m^3)
L = vol*10^(3); % liter of a Yeast(L)
K= 2*10^8*1000; % carrying capacity for 1L

% Simulation

sim = zeros(length(timev),1); % vector to store the values of N during simulation
i=0; % counter
max_t=0; % hours to fill 1L

for t=timev
    i=i+1;

    % main equation
    N = (K*N0*2^(t/tau2))/(K-N0+N0*2^(t/tau2));

    % find out how many hours it takes to fill 1L
    if max_t==0 if N>=K max_t=t; end; end

    % store value of N for plotting
    sim(i)=N;
end

% plot
figure(2); clf
plot(timev,sim,'-o')
xlabel('time (hours)');ylabel('population (cells)')
title('Logistic equation model for Yeast')

















