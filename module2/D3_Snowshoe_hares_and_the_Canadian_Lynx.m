% Simulation of Snowshoe hares and the Canadian Lynx population 
% considering in the predator-prey perspective

% The populations of both species have been observed to oscillate 
% with a cycle of about 10 years
% If under optimal conditions female hares bear an average of 18 young a year,
% 33% of which survive the rst month, the monthly survival rate of all hares 
% more than a month old is 95%, and sexual maturity is reached at one year

% AMS 333 Mathematical biology project 3
% David Hwang, 10/14/2024

%%%%%%%%%%%%%%%%%%%%%%%%
%% rate parameters of a hare(denoted as 'U') 
%%%%%%%%%%%%%%%%%%%%%%%%

Birth_U=(18*0.33*0.95^11)/2; % annual per capita reproduction rate (the young can survive to adult)
% Since males cannot produce offspring, we divide 2 here 
Death_U=1-1*0.95^12; % annual death rate of an originally existed hare(a parent)
Growth_U=1+Birth_U-Death_U; % a net annual per capita growth rate
alpha=log(Growth_U); % ln(R_a)

%%%%%%%%%%%%%%%%%%%%%%%%
%% rate parameters of Lynx(denoted as 'V')without hares
%%%%%%%%%%%%%%%%%%%%%%%%

Death_V=1-0.7^12; % annual death rate of an originally existing lynx
beta=-log(1-Death_V); % -ln(Ra_v)
% 1year=365days
gamma=1*365/1000;% annual per capita predation rate;assume that there are 1000 hares per square km^2
epsilon=1.5*0.1/10; % 10% of prey mass goes towards reproduction&rearing of kittens


%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation of Lokta-Volterra model 
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

U0=400; %initial value; 400 hares
V0=1; %initial value; 1 lynx
U=U0; % U is a number
V=V0; % V is a number

Birth_U=(18*0.33*0.95^11)/2; % annual per capita reproduction rate (the young can survive to adult)
Death_U=1-1*0.95^12; % annual death rate of an originally existed hare(a parent)
Growth_U=Birth_U+(1-Death_U); % a net annual per capita growth rate
alpha=log(Growth_U); % ln(Ra_u)

Death_V=1-0.7^12; % annual death rate of an originally existing lynx
beta=-log(1-Death_V); % -ln(Ra_v)

gamma=1*365/1000; % annual per capita predation rate;assume that there are 1000 hares per square km
epsilon=1.5*0.1/10; % 10% of prey mass goes towards reproduction&rearing of kittens


% simulation

dt=0.001; % time step; 0.001years
timev=0:dt:40; % time vector for 40 years which is 4 cycles
sim_u=zeros(length(timev),1); % vector to store the population of hares
sim_v=zeros(length(timev),1); % vector to store the population of lynx
k=0; %counter
sim_u(1)=U;
sim_v(1)=V;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        U = U + (alpha*U-gamma*U*V)*dt;
        V = V + (epsilon*gamma*U*V-beta*V)*dt;
    end

    if U<0 U=0; end;
    if V<0 V=0; end;
    
    % store value of U&V for plotting
    sim_u(k)=U;
    sim_v(k)=V;
end

% plot
figure(1); clf
hold on
plot(timev,sim_u,'LineWidth', 2)
plot(timev,sim_v,'LineWidth', 2)
hold off
legend('hare population','lynx population')
xlabel('Time (years)');ylabel('Population (numbers)')
title('Lokta-Volterra model')

figure(2); clf
hold on
plot(sim_u,sim_v,'LineWidth', 2)
plot([0,0],ylim,'-b','LineWidth', 1)
plot(0, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25); 
plot(xlim,[alpha/gamma,alpha/gamma],'-b','LineWidth', 1)
plot(xlim,[0,0],'-b','LineWidth', 1)
plot([beta/(epsilon*gamma),beta/(epsilon*gamma)],ylim,'-b','LineWidth', 1)
plot(beta/(epsilon*gamma), (alpha/gamma), 'r.', 'LineWidth', 2, 'MarkerSize', 25);
legend('population path','nullclines','stationary points')
xlabel('Number of hare (prey)');ylabel('Number of lynx (predator)')

hold off
title('Lynx-Hare phase plane')


%%%%%%%%%%%%%%%%%%%%%%%%
%% Observation 
%%%%%%%%%%%%%%%%%%%%%%%%

% amplitude
amp_u=(max(sim_u)-min(sim_u))/2; % amplitude of hare population
amp_v=(max(sim_v)-min(sim_v))/2; % amplitude of lynx population


% period of cycle of LynxHare phase plane

m_u=[]; % list for time of local maximum point of hare population
m_v=[]; % list for time of local maximum point of lynx population

k=0; % counter
for t=timev
    k=k+1;
    if sim_u(k)<sim_u(k+1)&&sim_u(k+2)<sim_u(k+1) m_u(end+1)=t; end;
    if sim_v(k)<sim_v(k+1)&&sim_v(k+2)<sim_v(k+1) m_v(end+1)=t; end;
    if k+2==length(sim_v) break; end;
end

c_u=m_u(2)-m_u(1); % cycle of hare pop
c_v=m_v(2)-m_v(1); % cycle of lynx pop
logic = c_u==c_v; % Through the phase plane and this, we can know the the population of lynx and the population of hare have the same cycle length.

% maximum and minimum population density
max_u=max(sim_u); % maximum population density of hare
min_u=min(sim_u); % minimum popuation density of hare
max_v=max(sim_v); % maximum population density of lynx
min_v=min(sim_v); % minimum population density of lynx

% maximum rate of change of population density
mrate_u=max(alpha.*sim_u-gamma.*sim_u.*sim_v); % maximum rate of hare
mrate_v=max(epsilon*gamma.*sim_u.*sim_v-beta.*sim_v); % maximum rate of lynx

% relative timing of the peak population for each species
relative_t = m_u-m_v; % always hare's population increase first and then lynx's population follow up and since the period of cycle is the same, the time difference of the peak time is the same throughout the whole simulation

% length of time that any species spends at below 1 individual per 100 km2

c1=0; % counter of how many times that lynx spends at below 1 individual per 100 km2 in one cycle.
k=m_v(1)/0.001; % counter starting from the first peak of hare's population
for t=m_v(1):dt:m_v(2)
    k=k+1;
    if 100*sim_v(k)<1 c1=c1+1; end;
end

c2 = c1*0.001; % the real time length of lynx spends at below 1 individual per 100 km2 in one cycle.
% min(sim_v)*100 is bigger than 1. so there is no time that lynx population
% below than 1.

% Since even in the 1km^2, the population of hare never goes down under 300
% population, we can know that the length of time that hare spens at below 1
% individual per 100km^2 is zero


%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation with different initial value 
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

U0=200; %initial value; 
V0=0.5; %initial value; 
U=U0; % U is a number
V=V0; % V is a number

Birth_U=(18*0.33*0.95^11)/2; % annual per capita reproduction rate (the young can survive to adult)
Death_U=1-1*0.95^12; % annual death rate of an originally existed hare(a parent)
Growth_U=Birth_U+(1-Death_U); % a net annual per capita growth rate
alpha=log(Growth_U); % ln(Ra_u)

Death_V=1-0.7^12; % annual death rate of an originally existing lynx
beta=-log(1-Death_V); % -ln(Ra_v)

gamma=1*365/1000; % annual per capita predation rate;assume that there are 1000 hares per square km
epsilon=1.5*0.1/10; % 10% of prey mass goes towards reproduction&rearing of kittens


% simulation

dt=0.001; % time step; 0.001years
timev=0:dt:40; % time vector for 40 years which is 4 cycles
sim_u=zeros(length(timev),1); % vector to store the population of hares
sim_v=zeros(length(timev),1); % vector to store the population of lynx
k=0; %counter
sim_u(1)=U;
sim_v(1)=V;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        U = U + (alpha*U-gamma*U*V)*dt;
        V = V + (epsilon*gamma*U*V-beta*V)*dt;
    end

    if U<0 U=0; end;
    if V<0 V=0; end;
    
    % store value of U&V for plotting
    sim_u(k)=U;
    sim_v(k)=V;
end

% plot
figure(1); clf
hold on
plot(timev,sim_u,'LineWidth', 2)
plot(timev,sim_v,'LineWidth', 2)
hold off
legend('hare population','lynx population')
xlabel('Time (years)');ylabel('Population (numbers)')
title('Lokta-Volterra model starts with prey 200 predator 0.5')

figure(2); clf
hold on
plot(sim_u,sim_v,'LineWidth', 2)
xline(0,'-b','LineWidth', 1)
plot(0, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25); 
yline(alpha/gamma,'-b','LineWidth', 1)
yline(0,'-b','LineWidth', 1)
xline(beta/(epsilon*gamma),'-b','LineWidth', 1)
plot(beta/(epsilon*gamma), (alpha/gamma), 'r.', 'LineWidth', 2, 'MarkerSize', 25);
legend('population path','nullclines','stationary points')
xlabel('Number of hare (prey)');ylabel('Number of lynx (predator)')
hold off
title('Lynx-Hare phase plane starts with prey 200 predator 0.5')



% when the population is 800/2
% => Since the population of lynx starts with big value, the hare population is not fluctuating much compared to when 400/1 
% So the size of phase plane? is small

% when the population is 200/.5
% => Since the population of lynx starts with small value, the hare
% population is fluctuating a lot.
% so the size of phase plane? is large


%%%%%%%%%%%%%%%%%%%%%%%%
%%  Lotka-Volterra model with logistic growth 
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

U0=400; %initial value; 
V0=1; %initial value; 
U=U0; % U is a number
V=V0; % V is a number

Birth_U=(18*0.33*0.95^11)/2; % annual per capita reproduction rate (the young can survive to adult)
Death_U=1-1*0.95^12; % annual death rate of an originally existed hare(a parent)
Growth_U=Birth_U+(1-Death_U); % a net annual per capita growth rate
alpha=log(Growth_U); % ln(Ra_u)
K=3000; % carrying capacity of 3000 per km^2
Death_V=1-0.7^12; % annual death rate of an originally existing lynx
beta=-log(1-Death_V); % -ln(Ra_v)

gamma=1*365/1000; % annual per capita predation rate;assume that there are 1000 hares per square km
epsilon=1.5*0.1/10; % 10% of prey mass goes towards reproduction&rearing of kittens


% simulation

dt=0.001; % time step; 0.001years
timev=0:dt:40; % time vector for 40 years which is 4 cycles
sim_u=zeros(length(timev),1); % vector to store the population of hares
sim_v=zeros(length(timev),1); % vector to store the population of lynx
i=0; %counter
sim_u(1)=U;
sim_v(1)=V;

for t=timev
    i=i+1;

    % main equation
    if i~=1
        U = U + (alpha*U*(1-U/K)-gamma*U*V)*dt;
        V = V + (epsilon*gamma*U*V-beta*V)*dt;
    end

    if U<0 U=0; end;
    if V<0 V=0; end;
    
    % store value of U&V for plotting
    sim_u(i)=U;
    sim_v(i)=V;
end

% plot
figure(1); clf
hold on
plot(timev,sim_u,'LineWidth', 2)
plot(timev,sim_v,'LineWidth', 2)
hold off
legend('hare population','lynx population')
xlabel('Time (years)');ylabel('Population (numbers)')
title('Lotka-Volterra model with logistic growth')

figure(2); clf
hold on
plot(sim_u,sim_v,'LineWidth', 2)
xline(0,'-b','LineWidth', 1)
plot(0, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25); 
plot(xlim,alpha/gamma*(1-xlim./K),'-b','LineWidth', 1)
yline(0,'-b','LineWidth', 1)
xline(beta/(epsilon*gamma),'-b','LineWidth', 1)
plot(beta/(epsilon*gamma), alpha/gamma*(1-beta/(epsilon*gamma*K)), 'r.', 'LineWidth', 2, 'MarkerSize', 25);
legend('population path','nullclines','stationary points')
xlabel('Number of hare (prey)');ylabel('Number of lynx (predator)')
hold off
title('Lynx-Hare phase plane with logistic growth')

% Compared to the previous original Lotka-Volterra model, we can notice
% from the phase plane and the pop VS. time graph that each populations are
% converges toward a stationary point(which is a stable stationary point)
% Similarities: each populations fluctuate at the beginning of the
% simulation
% Differences: As time goes by, we can find that this model goes to a
% stable stationary point. not-like a previous model which is just
% oscilating

%%%%%%%%%%%%%%%%%%%%%%%%
%% Lotka-Volterra model with Holling's disk equation
%%%%%%%%%%%%%%%%%%%%%%%%


% main parameters

U0=400; %initial value; 400 hares
V0=1; % initial value; 1 lynx
U=U0; % U is a number
V=V0; % V is a number

Birth_U=(18*0.33*0.95^11)/2; % annual per capita reproduction rate (the young can survive to adult)
Death_U=1-1*0.95^12; % annual death rate of an originally existed hare(a parent)
Growth_U=Birth_U+(1-Death_U); % a net annual per capita growth rate
alpha=log(Growth_U); % ln(Ra_u)

Death_V=1-0.7^12; % annual death rate of an originally existing lynx
beta=-log(1-Death_V); % -ln(Ra_v)

gamma=1*365/1000; % annual per capita predation rate;assume that there are 1000 hares per square km
epsilon=1.5*0.1/10; % 10% of prey mass goes towards reproduction&rearing of kittens
kappa = 4/24/365; % pre handling time of 4 hours

% simulation

dt=0.001; % time step; 0.001 years
timev=0:dt:40; % time vector for 40 years which is 4 cycles
sim_u=zeros(length(timev),1); % vector to store the population of hares
sim_v=zeros(length(timev),1); % vector to store the population of lynx
k=0; %counter
sim_u(1)=U;
sim_v(1)=V;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        U = U + (alpha*U-((gamma*U)/(1+gamma*kappa*U))*V)*dt;
        V = V + (epsilon*((gamma*U)/(1+gamma*kappa*U))*V-beta*V)*dt;
    end

    if U<0 U=0; end;
    if V<0 V=0; end;
    
    % store value of U&V for plotting
    sim_u(k)=U;
    sim_v(k)=V;
end

% plot
figure(1); clf
hold on
plot(timev,sim_u,'LineWidth', 2)
plot(timev,sim_v,'LineWidth', 2)
hold off
legend('hare population','lynx population')
xlabel('Time (years)');ylabel('Population (numbers)')
title("Lotka-Volterra model with Holling's disk equation")


x=0:1:9000; % x vector for plotting on the phase plane


figure(2); clf
hold on
plot(sim_u,sim_v,'LineWidth', 2)
xline(0,'-b','LineWidth', 1)
plot(0, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25); 
plot(x,alpha/gamma*(1+gamma*kappa.*x),'-b','LineWidth', 1)
yline(0,'-b','LineWidth', 1)
xline(beta/((epsilon-beta*kappa)*gamma),'-b','LineWidth', 1)
plot(beta/((epsilon-beta*kappa)*gamma), alpha/gamma*(1+gamma*kappa*(beta/((epsilon-beta*kappa)*gamma))), 'r.', 'LineWidth', 2, 'MarkerSize', 25);
legend('population path','nullclines','stationary points')
xlabel('Number of hare (prey)');ylabel('Number of lynx (predator)')
hold off
title('Final modified Lynx-Hare phase plane')






% since the hare's growth is sustained exponentially, if we modified the
% model considering the consuming limit of lynx which decreases the predataion rate, the amplitude of population of lynx
% keeps growing as time goes by.


%%%%%%%%%%%%%%%%%%%%%%%%
%% Final modified model with logistic and Holling's dist equation 
%%%%%%%%%%%%%%%%%%%%%%%%



% main parameters

U0=400; %initial value; 400 hares
V0=1; % initial value; 1 lynx
U=U0; % U is a number
V=V0; % V is a number

Birth_U=(18*0.33*0.95^11)/2; % annual per capita reproduction rate (the young can survive to adult)
Death_U=1-1*0.95^12; % annual death rate of an originally existed hare(a parent)
Growth_U=Birth_U+(1-Death_U); % a net annual per capita growth rate
alpha=log(Growth_U); % ln(Ra_u)
K=3000 ; % carrying capacity of 3000 per km^2

Death_V=1-0.7^12; % annual death rate of an originally existing lynx
beta=-log(1-Death_V); % -ln(Ra_v)

gamma=1*365/1000; % annual per capita predation rate;assume that there are 1000 hares per square km
epsilon=1.5*0.1/10; % 10% of prey mass goes towards reproduction&rearing of kittens
kappa = 4/(24*365); % pre handling time of 4 hours

% simulation

dt=0.001; % time step; 0.001years
timev=0:dt:40; % time vector for 40 years which is 4 cycles
sim_u=zeros(length(timev),1); % vector to store the population of hares
sim_v=zeros(length(timev),1); % vector to store the population of lynx
k=0; %counter
sim_u(1)=U;
sim_v(1)=V;

for t=timev
    k=k+1;

    % main equation
    if k~=1
        U = U + (alpha*U*(1-U/K)-((gamma*U)/(1+gamma*kappa*U))*V)*dt;
        V = V + (epsilon*((gamma*U)/(1+gamma*kappa*U))*V-beta*V)*dt;
    end

    if U<0 U=0; end;
    if V<0 V=0; end;
    
    % store value of U&V for plotting
    sim_u(k)=U;
    sim_v(k)=V;
end

x=0:1:1600; % x vector for plotting on the phase plane


% plot
figure(1); clf
hold on
plot(timev,sim_u,'LineWidth', 2)
plot(timev,sim_v,'LineWidth', 2)
hold off
legend('hare population','lynx population')
xlabel('Time (years)');ylabel('Population (numbers)')
title('Final modified Lokta-Volterra model')

figure(2); clf
hold on
plot(sim_u,sim_v,'LineWidth', 2)
xline(0,'-b','LineWidth', 1)
plot(0, 0, 'r.', 'LineWidth', 2, 'MarkerSize', 25); 
plot(x,alpha/gamma*(1-((1-gamma*kappa*K)/K).*x-(gamma*kappa/K).*x.^2),'-b','LineWidth', 1)
yline(0,'-b','LineWidth', 1)
xline(beta/((epsilon-beta*kappa)*gamma),'-b','LineWidth', 1)
plot(beta/((epsilon-beta*kappa)*gamma), alpha/gamma*(1-(1-gamma*kappa*K)/K*(beta/((epsilon-beta*kappa)*gamma))-(gamma*kappa/K)*(beta/((epsilon-beta*kappa)*gamma))^2), 'r.', 'LineWidth', 2, 'MarkerSize', 25);
legend('population path','nullclines','stationary points')
xlabel('Number of hare (prey)');ylabel('Number of lynx (predator)')
hold off
title('Final modified Lynx-Hare phase plane')



% In the realistic condition, there is a limited prey growth rate because of the limited resource and limited space etc.  
% And also there is a constraints to the predatation rate because of the
% limit of predator's consuming. 




%%


% period of cycle of LynxHare phase plane

m_u=[]; % list for time of local maximum point of hare population
m_v=[]; % list for time of local maximum point of lynx population

k=0; % counter
for t=timev
    k=k+1;
    if sim_u(k)<sim_u(k+1)&&sim_u(k+2)<sim_u(k+1) m_u(end+1)=t; end;
    if sim_v(k)<sim_v(k+1)&&sim_v(k+2)<sim_v(k+1) m_v(end+1)=t; end;
    if k+2==length(sim_v) break; end;
end

c_u=m_u(2)-m_u(1); % cycle of hare pop
c_v=m_v(2)-m_v(1); % cycle of lynx pop
logic = c_u==c_v; % Through the phase plane and this, we can know the the population of lynx and the population of hare have the same cycle length.



