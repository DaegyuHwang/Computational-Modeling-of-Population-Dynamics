% Simulation of malaria transmission model
% considering in the relation of Human and vector-borne, mosquito, perspective


% Compared to the cycle of malaria, houman life-cycle is long, we simplify
% our model by ignoring demographic processes amoung people. Bue we include
% birth and death rate of mosquitoes.
% We instroduce SIS model for human and SI model for mosquito.

% AMS 333 Mathematical biology project 4
% David Hwang, 12/12/2024


%%%%%%%%%%%%%%%%%%%%%%%%
%% SIS model for malaria (begin with no infective population for both species)
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

alpha = 0.02; % mosquito birth/death rate
beta = 5; % mosquito feeding rate
epsilon1 = 0.01; % fraction of bites that lead to infection(mosquito to human)
epsilon2 = 0.1; % fraction of bites that lead to infection(human to mosquito)
gamma = 0.05; % human recovery rate
total_Hp = 10^6; % total population of human
total_Mp = 10^6; % total population of mosquito
I_I_h = 0; % initial infective human population
I_S_h = total_Hp-I_I_h; % initial susceptive human population
I_I_m = 0; % initial infective mosquito population
I_S_m = total_Mp-I_I_m; % initial susceptive mosquito population


% simulation

dt=0.01; % time step: (one day/100)
timev=0:dt:200; % time vector for 200 days
sim_S_h=zeros(length(timev),1); % vector to store susceptive human population
sim_I_h=zeros(length(timev),1); % vector to store infective human population
sim_S_m=zeros(length(timev),1); % vector to store susceptive mosquitoe population
sim_I_m=zeros(length(timev),1); % vector to store infective mosquitoe population

k=0; %counter
S_H=I_S_h; % a number
I_H=I_I_h; % a number
S_M=I_S_m; % a number
I_M=I_I_m; % a number

for t=timev
    k=k+1;

    % main equation
    if k~=1
        S_H = S_H + (-epsilon1*beta*S_H/total_Hp*I_M + gamma*I_H)*dt;
        I_H = I_H + (epsilon1*beta*S_H/total_Hp*I_M - gamma*I_H)*dt;
        S_M = S_M + (alpha*total_Mp - epsilon2*beta*I_H/total_Hp*S_M - alpha*S_M)*dt;    
        I_M = I_M + (epsilon2*beta*I_H/total_Hp*S_M-alpha*I_M)*dt;
    end

    if S_H<0 S_H=0; end; % population cannot go below than zero
    if I_H<0 I_H=0; end;
    if S_M<0 S_M=0; end;
    if S_M<0 S_M=0; end;
    
    % store values for plotting
    sim_S_h(k)=S_H;
    sim_I_h(k)=I_H;
    sim_S_m(k)=S_M;
    sim_I_m(k)=I_M;
end

% plot
figure(1); clf
sgtitle('Starting from zero infective population for both species')
subplot(2,1,1)
hold on
plot(timev,sim_S_h,'LineWidth', 2)
plot(timev,sim_I_h,'LineWidth', 2)
hold off
legend('Susceptible(Human)','Infective(Human) ')
xlabel('Time (days)');ylabel('Population (capita)')
title('SIS model of Human')

subplot(2,1,2)
hold on
plot(timev,sim_S_m,'LineWidth', 2)
plot(timev,sim_I_m,'LineWidth', 2)
hold off
legend('Susceptible(Mosquito)','Infective(Mosquito)')
xlabel('Time (days)');ylabel('Population (capita)')
title('SI model of Mosquito')





%%%%%%%%%%%%%%%%%%%%%%%%
%% SIS model for malaria (only with initial infective mosquitoes)
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

alpha = 0.02; % mosquito birth/death rate
beta = 5; % mosquito feeding rate
epsilon1 = 0.01; % fraction of bites that lead to infection(mosquito to human)
epsilon2 = 0.1; % fraction of bites that lead to infection(human to mosquito)
gamma = 0.05; % human recovery rate
total_Hp = 10^6; % total population of human
total_Mp = 10^6; % total population of mosquito
I_I_h = 0; % initial infective human population
I_S_h = total_Hp-I_I_h; % initial susceptive human population
I_I_m = 0; % initial infective mosquito population
I_S_m = total_Mp-I_I_m; % initial susceptive mosquito population


% simulation

dt=0.01; % time step: (one day/100)
timev=0:dt:200; % time vector for 200 days
sim_S_h=zeros(length(timev),1); % vector to store susceptive human population
sim_I_h=zeros(length(timev),1); % vector to store infective human population
sim_S_m=zeros(length(timev),1); % vector to store susceptive mosquitoe population
sim_I_m=zeros(length(timev),1); % vector to store infective mosquitoe population

k=0; j=0; l=-1;%counter
S_H=I_S_h; % a number
I_H=I_I_h; % a number
S_M=I_S_m; % a number
I_M=I_I_m; % a number




figure(2); clf;
sgtitle('Started with only infected mosquitoes')

for i = 0.1:0.4:0.9
    j=j+2;l=l+2;k=0;
    S_H=I_S_h; % a number
    I_H=I_I_h; % a number
    I_M = i*total_Mp;
    S_M = total_Mp - I_M;

    for t=timev
        k=k+1;
    
        % main equation
        if k~=1
            S_H = S_H + (-epsilon1*beta*S_H/total_Hp*I_M + gamma*I_H)*dt;
            I_H = I_H + (epsilon1*beta*S_H/total_Hp*I_M - gamma*I_H)*dt;
            S_M = S_M + (alpha*total_Mp - epsilon2*beta*I_H/total_Hp*S_M - alpha*S_M)*dt;    
            I_M = I_M + (epsilon2*beta*I_H/total_Hp*S_M-alpha*I_M)*dt;
        end
    
        if S_H<0 S_H=0; end; % population cannot go below than zero
        if I_H<0 I_H=0; end;
        if S_M<0 S_M=0; end;
        if S_M<0 S_M=0; end;
        
        % store values for plotting
        sim_S_h(k)=S_H;
        sim_I_h(k)=I_H;
        sim_S_m(k)=S_M;
        sim_I_m(k)=I_M;
    end

    % plot 
    
    
    subplot(3,2,l)
    hold on
    plot(timev,sim_S_h(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_h(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Human)','Infective(Human) ')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SIS model of Human (with %f of infective mosquito)', i))
    
    subplot(3,2,j)
    hold on
    plot(timev,sim_S_m(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_m(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Mosquito)','Infective(Mosquito)')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SI model of Mosquito (with %f of infective mosquito)', i))
    
         
end

 

%%%%%%%%%%%%%%%%%%%%%%%%
%% SIS model for malaria (only with initial infectived human population)
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

alpha = 0.02; % mosquito birth/death rate
beta = 5; % mosquito feeding rate
epsilon1 = 0.01; % fraction of bites that lead to infection(mosquito to human)
epsilon2 = 0.1; % fraction of bites that lead to infection(human to mosquito)
gamma = 0.05; % human recovery rate
total_Hp = 10^6; % total population of human
total_Mp = 10^6; % total population of mosquito
I_I_h = 0; % initial infective human population
I_S_h = total_Hp-I_I_h; % initial susceptive human population
I_I_m = 0; % initial infective mosquito population
I_S_m = total_Mp-I_I_m; % initial susceptive mosquito population


% simulation

dt=0.01; % time step: (one day/100)
timev=0:dt:200; % time vector for 200 days
sim_S_h=zeros(length(timev),1); % vector to store susceptive human population
sim_I_h=zeros(length(timev),1); % vector to store infective human population
sim_S_m=zeros(length(timev),1); % vector to store susceptive mosquitoe population
sim_I_m=zeros(length(timev),1); % vector to store infective mosquitoe population

k=0; j=0; l=-1;%counter
S_H=I_S_h; % a number
I_H=I_I_h; % a number
S_M=I_S_m; % a number
I_M=I_I_m; % a number




figure(3); clf;
sgtitle('Started with only infected human')

for i = 0.1:0.4:0.9
    j=j+2;l=l+2;k=0;

    I_H=i*total_Hp; 
    S_H = total_Hp - I_H;
    S_M=I_S_m; % a number
    I_M=I_I_m; % a number

    for t=timev
        k=k+1;
    
        % main equation
        if k~=1
            S_H = S_H + (-epsilon1*beta*S_H/total_Hp*I_M + gamma*I_H)*dt;
            I_H = I_H + (epsilon1*beta*S_H/total_Hp*I_M - gamma*I_H)*dt;
            S_M = S_M + (alpha*total_Mp - epsilon2*beta*I_H/total_Hp*S_M - alpha*S_M)*dt;    
            I_M = I_M + (epsilon2*beta*I_H/total_Hp*S_M-alpha*I_M)*dt;
        end
    
        if S_H<0 S_H=0; end; % population cannot go below than zero
        if I_H<0 I_H=0; end;
        if S_M<0 S_M=0; end;
        if S_M<0 S_M=0; end;
        
        % store values for plotting
        sim_S_h(k)=S_H;
        sim_I_h(k)=I_H;
        sim_S_m(k)=S_M;
        sim_I_m(k)=I_M;
    end

    % plot 
    
    
    subplot(3,2,l)
    hold on
    plot(timev,sim_S_h(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_h(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Human)','Infective(Human) ')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SIS model of Human (with %f of infective human)', i))
    
    subplot(3,2,j)
    hold on
    plot(timev,sim_S_m(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_m(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Mosquito)','Infective(Mosquito)')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SI model of Mosquito (with %f of infective human)', i))
    
         
end

 
%%%%%%%%%%%%%%%%%%%%%%%%
%% SIS model for malaria (run with different beta values)
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

alpha = 0.02; % mosquito birth/death rate
beta_list = [5,2,1,0.5]; % mosquito feeding rate list
epsilon1 = 0.01; % fraction of bites that lead to infection(mosquito to human)
epsilon2 = 0.1; % fraction of bites that lead to infection(human to mosquito)
gamma = 0.05; % human recovery rate
total_Hp = 10^6; % total population of human
total_Mp = 10^6; % total population of mosquito
I_I_h = 0.1*total_Hp; % initial infective human population(10%)
I_S_h = total_Hp-I_I_h; % initial susceptive human population
I_I_m = 0.1*total_Mp; % initial infective mosquito population(10%)
I_S_m = total_Mp-I_I_m; % initial susceptive mosquito population


% simulation

dt=0.01; % time step: (one day/100)
timev=0:dt:200; % time vector for 200 days
sim_S_h=zeros(length(timev),1); % vector to store susceptive human population
sim_I_h=zeros(length(timev),1); % vector to store infective human population
sim_S_m=zeros(length(timev),1); % vector to store susceptive mosquitoe population
sim_I_m=zeros(length(timev),1); % vector to store infective mosquitoe population

k=0; j=0; l=-1;%counter
S_H=I_S_h; % a number
I_H=I_I_h; % a number
S_M=I_S_m; % a number
I_M=I_I_m; % a number




figure(4); clf;
sgtitle('Simulation with different beta values')

for i = 1:1:4
    j=j+2;l=l+2;k=0;
    beta = beta_list(i);

    S_H=I_S_h; % a number
    I_H=I_I_h; % a number
    S_M=I_S_m; % a number
    I_M=I_I_m; % a number

    for t=timev
        k=k+1;
    
        % main equation
        if k~=1
            S_H = S_H + (-epsilon1*beta*S_H/total_Hp*I_M + gamma*I_H)*dt;
            I_H = I_H + (epsilon1*beta*S_H/total_Hp*I_M - gamma*I_H)*dt;
            S_M = S_M + (alpha*total_Mp - epsilon2*beta*I_H/total_Hp*S_M - alpha*S_M)*dt;    
            I_M = I_M + (epsilon2*beta*I_H/total_Hp*S_M-alpha*I_M)*dt;
        end
    
        if S_H<0 S_H=0; end; % population cannot go below than zero
        if I_H<0 I_H=0; end;
        if S_M<0 S_M=0; end;
        if S_M<0 S_M=0; end;
        
        % store values for plotting
        sim_S_h(k)=S_H;
        sim_I_h(k)=I_H;
        sim_S_m(k)=S_M;
        sim_I_m(k)=I_M;
    end

    % plot 
    
    
    subplot(4,2,l)
    hold on
    plot(timev,sim_S_h(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_h(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Human)','Infective(Human) ')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SIS model of Human (with beta = %f)', beta))
    
    subplot(4,2,j)
    hold on
    plot(timev,sim_S_m(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_m(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Mosquito)','Infective(Mosquito)')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SI model of Mosquito (with beta = %f)', beta))
    
         
end


figure(10);clf;
plot(sim_I_m,sim_S_m)

%%%%%%%%%%%%%%%%%%%%%%%%
%% SIS model for malaria (run with different total population)
%%%%%%%%%%%%%%%%%%%%%%%%

% main parameters

alpha = 0.02; % mosquito birth/death rate
beta = 5; % mosquito feeding rate list
epsilon1 = 0.01; % fraction of bites that lead to infection(mosquito to human)
epsilon2 = 0.1; % fraction of bites that lead to infection(human to mosquito)
gamma = 0.05; % human recovery rate
total_Hp = 10^6; % total population of human
total_Mp_list = [10^6, 10^5, 10^4]; % total population of mosquito
I_I_h = 0.1*total_Hp; % initial infective human population(10%)
I_S_h = total_Hp-I_I_h; % initial susceptive human population
%I_I_m = 0.1*total_Mp; % initial infective mosquito population(10%)
%I_S_m = total_Mp-I_I_m; % initial susceptive mosquito population


% simulation

dt=0.01; % time step: (one day/100)
timev=0:dt:500; % time vector for 200 days
sim_S_h=zeros(length(timev),1); % vector to store susceptive human population
sim_I_h=zeros(length(timev),1); % vector to store infective human population
sim_S_m=zeros(length(timev),1); % vector to store susceptive mosquitoe population
sim_I_m=zeros(length(timev),1); % vector to store infective mosquitoe population

k=0; j=0; l=-1;%counter
S_H=I_S_h; % a number
I_H=I_I_h; % a number
S_M=I_S_m; % a number
I_M=I_I_m; % a number




figure(5); clf;
sgtitle('Simulation with different mosquito total population')

for i = 1:1:3
    j=j+2;l=l+2;k=0;
    total_Mp = total_Mp_list(i);
    I_I_m = 0.1*total_Mp;  
    I_S_m = total_Mp-I_I_m;  

    S_H=I_S_h; % a number
    I_H=I_I_h; % a number
    S_M=I_S_m; % a number
    I_M=I_I_m; % a number

    for t=timev
        k=k+1;
    
        % main equation
        if k~=1
            S_H = S_H + (-epsilon1*beta*S_H/total_Hp*I_M + gamma*I_H)*dt;
            I_H = I_H + (epsilon1*beta*S_H/total_Hp*I_M - gamma*I_H)*dt;
            S_M = S_M + (alpha*total_Mp - epsilon2*beta*I_H/total_Hp*S_M - alpha*S_M)*dt;    
            I_M = I_M + (epsilon2*beta*I_H/total_Hp*S_M-alpha*I_M)*dt;
        end
    
        if S_H<0 S_H=0; end; % population cannot go below than zero
        if I_H<0 I_H=0; end;
        if S_M<0 S_M=0; end;
        if S_M<0 S_M=0; end;
        
        % store values for plotting
        sim_S_h(k)=S_H;
        sim_I_h(k)=I_H;
        sim_S_m(k)=S_M;
        sim_I_m(k)=I_M;
    end

    % plot 
    
    
    subplot(3,2,l)
    hold on
    plot(timev,sim_S_h(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_h(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Human)','Infective(Human) ')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SIS model of Human (with mosqutio population = %f)', total_Mp))
    
    subplot(3,2,j)
    hold on
    plot(timev,sim_S_m(1:length(timev)),'LineWidth', 2)
    plot(timev,sim_I_m(1:length(timev)),'LineWidth', 2)
    hold off
    legend('Susceptible(Mosquito)','Infective(Mosquito)')
    xlabel('Time (days)');ylabel('Population (capita)')
    title(sprintf('SI model of Mosquito (with mosquito population = %f)', total_Mp))
    
         
end












