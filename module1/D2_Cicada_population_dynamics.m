% Simulation of Cicada population dynamics model (periodic growth cycle)

% If, in absence of competition, a cicada population is equally split between males and females,
% 80% of adults are able to reproduce, each female lays 300 eggs, and 20% of eggs survive to adulthood
% Since the life cycle of cicada is 17year, we set one cycle period is 17y

% David Hwang, 09/19/24

%%%%%%%%%%%%%%%%%%%%
%% finding carrying capacity and R0
%%%%%%%%%%%%%%%%%%%%

% the life cycle of cicada is 17y, so we assume that the populations in one cycle are
% not added to the next cycle (All of them only can survive for one cycle)

K1 = 250000000; % carrying capacity (km^2)
K = K1*350; % carrying capacity of 350km^2 region

R0 = (1/2)*0.8*300*0.2; % growth rate per capita(in terms of the total population)



%%%%%%%%%%%%%%%%%%%%
%% Discrete growth model(Hassell equation)
%%%%%%%%%%%%%%%%%%%%

% main parameters

N0 = 100; % initial population of cicada
N = N0; % N is a number; N0 is a initial value 
R0 = (1/2)*0.8*300*0.2; % growth rate per capita
K = 250000000*350; % carrying capacity of 350km^2 region
a = R0/K; % a parameter of Hassell equation
b = 1; % a parameter of Hassell equation

% simulation

timev = 0:17:187; % time vector for 170years which is 11 cycles
sim = zeros(length(timev),1); % vector to store the values of N during simulation
k=0; % counter

for t=timev
    k = k+1;

    % main equation
    if k~=1 N=(R0*N)/((1+a*N)^b); end

    % store value of N for plotting
    sim(k)=N;
end


% plot

figure(1); clf
hold on
stairs(sim(1:k-1),sim(2:k),'-o')
plot(xlim,ylim,'-b')
hold off
xlabel('N(n)');ylabel('N(n+1)')
title('N(n+1) VS N(n)')



%%%%%%%%%%%%%%%%%%%%
%% Discrete growth model(Hassell equation)
%%%%%%%%%%%%%%%%%%%%


% main parameters

N0 = 1000000000; % initial population of cicada
N = N0; % N is a number; N0 is a initial value 
R0 = (1/2)*0.8*300*0.2; % growth rate per capita
K = 250000000*350; % carrying capacity of 350km^2 region
a = R0/K; % a parameter of Hassell equation
b = 1; % a parameter of Hassell equation

% simulation

timev = 0:17:187; % time vector for 170years which is 11 cycles
sim = zeros(length(timev),1); % vector to store the values of N during simulation
k=0; % counter

for t=timev
    k = k+1;

    % main equation
    if k~=1 N=(R0*N)/((1+a*N)^b); end

    % store value of N for plotting
    sim(k)=N;
end


% plot

figure(1); clf
hold on
plot(timev,sim,'-o')
yline(K,'--r');
xlabel('time (years)');ylabel('population ')
title('population VS. year (when N0=1000000000)')

%%%%%%%%%%%%%%%%%%%%
%% Simulation about the amount of popultion
%%%%%%%%%%%%%%%%%%%%

% main parameters % we set the initial population to 100

N0 = 100; % initial population of cicada
N = N0; % N is a number; N0 is a initial value 
R0 = (1/2)*0.8*300*0.2; % growth rate per capita
K = 250000000*350; % carrying capacity of 350km^2 region
a = R0/K; % a parameter of Hassell equation
b = 1; % a parameter of Hassell equation


% simulation

timev = 0:17:187; % time vector for 170years which is 11 cycles
sim = zeros(length(timev),1); % vector to store the values of N during simulation
k=0; % counter

for t=timev
    k = k+1;

    % main equation
    if k~=1 N=(R0*N)/((1+a*N)^b); end

    % store value of N for plotting
    sim(k)=N;
end

total_g = sim*2; % populaion vector in respect to grams
total_kg = total_g/1000; % populaion vector in respect to kgrams
total_p = fix(total_kg/70); % vector for the number of people with an equivalent mass


%%%%%%%%%%%%%%%%%%%%
%% simulation about the time to recover the number of individuals
%%%%%%%%%%%%%%%%%%%%





% main parameters % we set the initial population to 100

N0 = 100000; % initial population of cicada
N = N0; % N is a number; N0 is a initial value 
R0 = (1/2)*0.8*300*0.2; % growth rate per capita
K = 250000000*350; % carrying capacity of 350km^2 region
a = R0/K; % a parameter of Hassell equation
b = 1; % a parameter of Hassell equation


% simulation

timev = 0:17:187; % time vector for 170years which is 11 cycles
sim = zeros(length(timev),1); % vector to store the values of N during simulation
k=0; % counter

for t=timev
    k = k+1;

    % main equation
    if k~=1 N=(R0*N)/((1+a*N)^b); end

    % store value of N for plotting
    sim(k)=N;
end

total_g = sim*2; % populaion vector in respect to grams
total_kg = total_g/1000; % populaion vector in respect to kgrams
total_p = fix(total_kg/70); % vector for the number of people with an equivalent mass



timev(7) % 102 years
total_kg(7) % total kg after 102 years
total_p(7) % the number of people after 102 years


