
% For COVID-19, estimates early in the pandemic suggested a basic reproductive
% factor of about 3.0, with an average length of infectivity of about 10 days.


%%%%%%%%%%%%%%%%%%%%%%%%
%% SIR model of COVID-19 
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

r_0 = 3; % basic reproductive factor
tau = 10; % average length of infectivity (days)
gamma = 1/tau; % rate of recovery
beta = r_0/tau; % likelihood of new infection
inf_0 = 4; % 4 infected individuals were arrived in a city
total_pop = 10000000+inf_0; % total number of population in a city

R_hat0 = 0; % initial fraction of removed individuals
I_hat0 = inf_0/total_pop; % initial fraction of infective individuals
S_hat0 = 1-I_hat0-R_hat0; % initial fraction of susceptible individuals


I_hat = I_hat0; % I_hat is a number
S_hat = S_hat0; % S_hat is a number
R_hat = R_hat0; % R_hat is a number


% simulation

dt=1; % time step: one day
timev=0:dt:300; % time vector for 40 years
sim_S_h=zeros(length(timev),1); % vector to store the fraction of susceptible
sim_I_h=zeros(length(timev),1); % vector to store the fraction of infective
sim_R_h=zeros(length(timev),1); % vector to store the fraction of removed

k=0; %counter
sim_S_h(1)=I_hat;
sim_I_h(1)=I_hat;
sim_R_h(1)=R_hat;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        S_hat = sim_S_h(k-1) - beta*sim_S_h(k-1)*sim_I_h(k-1);
        I_hat = I_hat + beta*sim_S_h(k-1)*sim_I_h(k-1) - gamma*sim_I_h(k-1);
        R_hat = R_hat + gamma*sim_I_h(k-1);         
    end

    if S_hat<0 S_hat=0; end;
    if R_hat>1 R_hat=1; end;
    
    % store value of S&I&R for plotting
    sim_S_h(k)=S_hat;
    sim_I_h(k)=I_hat;
    sim_R_h(k)=R_hat;

end

figure(1); clf
hold on
plot(timev,sim_S_h,'LineWidth', 2)
plot(timev,sim_I_h,'LineWidth', 2)
plot(timev,sim_R_h,'LineWidth', 2)
hold off
legend('Susceptible','Infective','Removed')
xlabel('Time (days)');ylabel('Fraction (target/total)')
title('SIR model of COVIE-19')

%%%%%%%%%%%%%%%%%%%%%%%%
%% Questions for insights
%%%%%%%%%%%%%%%%%%%%%%%%

% peak time of infected individuals
[max_I,idx_m] = max(sim_I_h)
% active cases per 100,000 people at the peak
max_I*100000

% 84 days after the arrival of the first infected individulas,
% the number of infected individuals becomes highest.
% Through multiply the peak fraction of infected individuals by 100,000
% , we can get the value for the number of active cases per
% 100,000 people, around 3.1235e+04.

%%
% the highest number of new daily infections
Daily_I = beta*sim_S_h.*sim_I_h; % daily inflow fraction of I
[max_d,idx_d] = max(Daily_I)
max_d*100000

% From the daily fraction or the infection occurence list(beta*S_hat*I_hat),
% we can get the maximum 0.0421 at the 78th index. Thus, the highest number of new daily
% infections occur between 77-78 days after the arrival of the first infected individulas.
% And this also means around 4.2108e+03 people get infected per 100,000 people at that day.

%%  
 
Daily_I = beta*sim_S_h.*sim_I_h; % daily inflow fraction of I
Daily_R = gamma*sim_I_h; % daily inflow fraction of R
dynmaic_I = Daily_I-Daily_R; % derivative of I
idx_epidemic = find(dynmaic_I<=0, 1, 'first') % index of the end of the epidemic
sim_I_h(idx_epidemic)*100000
(sim_I_h(idx_epidemic)+sim_R_h(idx_epidemic))*100000
 

%% 
% Fraction of infected population by the end of the epidemic
sim_I_h(end)
sim_R_h(end)
sim_S_h(end)
1-sim_S_h(end)

% we can find that the peidemic ends at 85th index which means after 84days
% the epidemic is finished and we can find that 3.1235e+04 people are
% currently infected and if we include 'R_hat' people, we can say that
% 6.8789e+04 people have been infected from the beginning to the end of
% this simulation.

%%
% need of hospital care at the peak of the epidemic assumed 10% of Hospitalization
max(sim_I_h)*0.1*100000 % needed of hospital beds at peak demand per 100,000people
max_d*0.1*100000 % expected hospital admissions at peak growth per 100,000people

% From the result, we can know that around 3.1235e+03 beds would be needed at
% peak demand and also when the growth peaked, we can expect that
% 421.0761 people would get hospital admission at that day.

%%
% # of people who died by the end of the epidemic assumed that the mortality of the illness was 1.5%
Death_rate = (I_hat0+sum(Daily_I))*0.015;
Death_rate*total_pop


% From the result we can know that 1.4200e+05 people would have died by the
% end of the epidemic. 

%%%%%%%%%%%%%%%%%%%%%%%%
%% Repeat simulation with different value of parameters
%%%%%%%%%%%%%%%%%%%%%%%%


% main parameters

r_0 = 3; % basic reproductive factor
tau = 10; % average length of infectivity (days)
gamma = 1/tau; % rate of recovery
beta = r_0/tau; % likelihood of new infection
inf_0 = 400; % 4 infected individuals were arrived in a city
total_pop = 10000000+inf_0; % total number of population in a city

R_hat0 = 0; % initial fraction of removed individuals
I_hat0 = inf_0/total_pop; % initial fraction of infective individuals
S_hat0 = 1-I_hat0-R_hat0; % initial fraction of susceptible individuals


I_hat = I_hat0; % I_hat is a number
S_hat = S_hat0; % S_hat is a number
R_hat = R_hat0; % R_hat is a number


% simulation

dt=1; % time step: one day
timev=0:dt:300; % time vector for 40 years
sim_S_h=zeros(length(timev),1); % vector to store the fraction of susceptible
sim_I_h=zeros(length(timev),1); % vector to store the fraction of infective
sim_R_h=zeros(length(timev),1); % vector to store the fraction of removed

k=0; %counter
sim_S_h(1)=I_hat;
sim_I_h(1)=I_hat;
sim_R_h(1)=R_hat;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        S_hat = sim_S_h(k-1) - beta*sim_S_h(k-1)*sim_I_h(k-1);
        I_hat = I_hat + beta*sim_S_h(k-1)*sim_I_h(k-1) - gamma*sim_I_h(k-1);
        R_hat = R_hat + gamma*sim_I_h(k-1);         
    end

    if S_hat<0 S_hat=0; end;
    if R_hat>1 R_hat=1; end;
    
    % store value of S&I&R for plotting
    sim_S_h(k)=S_hat;
    sim_I_h(k)=I_hat;
    sim_R_h(k)=R_hat;

end

figure(2); clf
hold on
plot(timev,sim_S_h,'LineWidth', 2)
plot(timev,sim_I_h,'LineWidth', 2)
plot(timev,sim_R_h,'LineWidth', 2)
hold off
legend('Susceptible','Infective','Removed')
xlabel('Time (days)');ylabel('Fraction (target/total)')
title('(iii) SIR model of COVIE-19 with I0=400')



%%%%%%%%%%%%%%%%%%%%%%%%
%% Repeat simulation with different value of r0
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

r_0 = 1; % basic reproductive factor
tau = 10; % average length of infectivity (days)
gamma = 1/tau; % rate of recovery
beta = r_0/tau; % likelihood of new infection
inf_0 = 4; % 4 infected individuals were arrived in a city
total_pop = 10000000+inf_0; % total number of population in a city

I_hat0 = inf_0/total_pop; % initial fraction of infective individuals
S_hat0 = 1-I_hat0; % initial fraction of susceptible individuals
R_hat0 = 0; % initial fraction of removed individuals

I_hat = I_hat0; % I_hat is a number
S_hat = S_hat0; % S_hat is a number
R_hat = R_hat0; % R_hat is a number


% simulation

dt=1; % time step: one day
timev=0:dt:600; % time vector for 40 years
sim_S_h=zeros(length(timev),1); % vector to store the fraction of susceptible
sim_I_h=zeros(length(timev),1); % vector to store the fraction of infective
sim_R_h=zeros(length(timev),1); % vector to store the fraction of removed

k=0; %counter
sim_S_h(1)=I_hat;
sim_I_h(1)=I_hat;
sim_R_h(1)=R_hat;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        S_hat = sim_S_h(k-1) - beta*sim_S_h(k-1)*sim_I_h(k-1);
        I_hat = I_hat + beta*sim_S_h(k-1)*sim_I_h(k-1) - gamma*sim_I_h(k-1);
        R_hat = R_hat + gamma*sim_I_h(k-1);         
    end

    if S_hat<0 S_hat=0; end;
    if R_hat>1 R_hat=1; end;
    
    % store value of S&I&R for plotting
    sim_S_h(k)=S_hat;
    sim_I_h(k)=I_hat;
    sim_R_h(k)=R_hat;

end

figure(3); clf
hold on
 
plot(timev,sim_I_h,'LineWidth', 2)
plot(timev,sim_R_h,'LineWidth', 2)
hold off
legend( 'Infective','Removed')
xlabel('Time (days)');ylabel('Fraction (target/total)')
title('(iii) SIR model of COVIE-19 with r0=0.5')

 

%%%%%%%%%%%%%%%%%%%%%%%%
%% finding R0 for "herd immunity"
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

r_0 = 3; % basic reproductive factor
tau = 10; % average length of infectivity (days)
gamma = 1/tau; % rate of recovery
beta = r_0/tau; % likelihood of new infection
inf_0 = 4; % 4 infected individuals were arrived in a city
total_pop = 10000000+inf_0; % total number of population in a city

R_hat0 = 0.75; % initial fraction of removed individuals
I_hat0 = inf_0/total_pop; % initial fraction of infective individuals
S_hat0 = 1-I_hat0-R_hat0; % initial fraction of susceptible individuals

I_hat = I_hat0; % I_hat is a number
S_hat = S_hat0; % S_hat is a number
R_hat = R_hat0; % R_hat is a number


% simulation

dt=1; % time step: one day
timev=0:dt:1000; % time vector for 40 years
sim_S_h=zeros(length(timev),1); % vector to store the fraction of susceptible
sim_I_h=zeros(length(timev),1); % vector to store the fraction of infective
sim_R_h=zeros(length(timev),1); % vector to store the fraction of removed

k=0; %counter
sim_S_h(1)=I_hat;
sim_I_h(1)=I_hat;
sim_R_h(1)=R_hat;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        S_hat = sim_S_h(k-1) - beta*sim_S_h(k-1)*sim_I_h(k-1);
        I_hat = I_hat + beta*sim_S_h(k-1)*sim_I_h(k-1) - gamma*sim_I_h(k-1);
        R_hat = R_hat + gamma*sim_I_h(k-1);         
    end

    if S_hat<0 S_hat=0; end;
    if R_hat>1 R_hat=1; end;
    
    % store value of S&I&R for plotting
    sim_S_h(k)=S_hat;
    sim_I_h(k)=I_hat;
    sim_R_h(k)=R_hat;

end

figure(4); clf
hold on
plot(timev,sim_S_h,'LineWidth', 2)
plot(timev,sim_I_h,'LineWidth', 2)
plot(timev,sim_R_h,'LineWidth', 2)
hold off
legend('Susceptible','Infective','Removed')
xlabel('Time (days)');ylabel('Fraction (target/total)')
title('(iii) SIR model of COVIE-19 with R0=0.75')

 
 