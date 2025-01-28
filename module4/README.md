#Malaria transmission model(Vector-borne disease)

This directory contains MATLAB code and the corresponding lab report that simulates and predicts the infectious population for the Malaria.

Purpose: Running mathematical models for the population dynamics of a disease, considering properties of the disease, to predict their growth of the disease in the real world and discuss the implications of the the results.

#Used Mathematical model

We introduced SIS (susceptible, infected, susceptible) model to simulate the impact of Malaria on humans, considering its property that individuals do not develop immunity after an infection. This means that individuals in the "infected" population return to the "susceptible" population when they recover from the disease. Most importantly, Malaria is one of the Vector-borne diseases so it can not be modeled by considering one species because they are not transmitted directly from individual to individual. Thus, we add another SI model to the system to describe the mosquito popultion (their lifespan is not long enough to return to the "susceptible" population after recovery). Additionally, we ignore demographic processes among people except the "births" and "deaths" of mosquitoes. Because we are interested in situations where a disease is introduced into a society where population dynamics have already converged to a point.

The Malaria transmission model (SIS & SI model):
![image](https://github.com/user-attachments/assets/2ca945e8-1b73-4b63-b1c2-66d8835463d7)

