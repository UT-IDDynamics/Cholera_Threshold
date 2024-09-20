# Quantifying the Impact of Expanded Antibiotic Treatment on Cholera Outbreaks: A Theoretical and Applied Framework

This project combines the theoretical and applied research from two studies to explore and quantify the effects of expanding antibiotic treatment to include moderately symptomatic cholera infections. 

## Overview

The project is divided into two main components:

1. **Quantifying the Impact of Expanded Antibiotic Treatment on Cholera Outbreaks:** This section focuses on assessing how expanding antibiotic treatment guidelines to include moderately symptomatic cholera infections could influence the overall burden of cholera and antibiotic stewardship. The primary goal is to quantify the impact of such a policy change on outbreak dynamics and antibiotic usage.

2. **A Theoretical Framework to Quantify the Tradeoff Between Individual and Population Benefits of Expanded Antibiotic Use:** This part of the project develops a mathematical framework to explore scenarios where expanded antibiotic treatment could be beneficial for both individual and population health. By analyzing the conditions under which treating moderate infections can decrease transmission and reduce the total number of antibiotic doses administered, we aim to provide a comprehensive understanding of the potential benefits and risks associated with this strategy.[Read the preprint here](https://www.medrxiv.org/content/10.1101/2024.08.28.24312731v1.full)

## Quantifying the Impact of Expanded Antibiotic Treatment on Cholera Outbreaks

by Sharia Ahmed, Cormac LaPrete, Iza Ciglenecki, Andrew Azman, Daniel Leung, Lindsay Keegan

To assess the impact of expanded antibiotic treatment guidelines to include moderately symptomatic infections on cholera outbreaks, we simulate cholera transmission in a non-endemic setting using a compartmental model (Figure 1). Our objective was to assess how such a policy change could impact the overall burden of cholera as well the impact on antibiotic stewardship. We evaluate the impact treating moderately symptomatic infections with antibiotics different outbreak scenarios (varying $$R_e$$), different treatment-seeking scenarios (varying the proportion of moderate infections who seek treatment), and different treatment guidelines (varying the proportion of treatment-seeking moderate infections who recieve antibiotics). 

<div align="center">
  <img src="images/CholeraEpiModel.png" alt="Figure1" width="500"/>
</div>

Figure 1: Compartmental diagram of the cholera transmission model. All individuals start as susceptible (S) and become exposed (E) at a rate $$\lambda$$. Exposed individuals transition to the infected (I) compartment and we differentiate by symptoms ($$I_A,I_M,I_S$$, asymptomatic, moderate or severely symptomatic, respectively) and by treatment seeking behavior $$(I_{U}, I_H)$$, not healthcare seeking and healthcare seeking, respectively). All severely symptomatic infections are treated with antibiotics whereas not all moderately symptomatic infections who seek treatment receive antibiotics. The proportion of healthcare seeking moderate infections who receive antibiotics is governed by q. Untreated infections, both moderate and severe continue to shed for a longer duration following the resolution of symptoms ($$I_{Msh},I_{Ssh}$$), occurring at rate $$\alpha_M,\alpha_S$$, respectively, moderately symptomatic infections who are treated with antibiotics continue to shed for a shorter duration following treatment ($$I_{Mabx}$$), occurring at rate $$\delta\theta$$, whereas severely symptomatic infections who are treated with antibiotics remain in a treatment facility and do not contribute to transmission ($$I_{Sabx}$$). Infectious individuals either recover (R, at rates $$\gamma_A,\gamma_M,\gamma_{Mabx},\gamma_S,\gamma_{Sabx}$$) or die (D, at rates $$\mu_M,\mu_S$$) and we differentiate between individuals who recover without antibiotic treatment ($$R_{un}$$) and those who recover with antibiotics ($$R_{abx}$$) to compare the number of doses used under different treatment scenarios.

Through simulation, we show that expanding treatment guidelines to include moderately symptomatic infections can substantially reduce the burden of cholera in low and intermediate transmission settings, especially when rates of treatment-seeking behavior is high. Further, in some scenarios, we find that expanding antibiotic treatment criteria provides a public health benefit, as each additional dose averts more than one infection and in limited scenarios, we showed that expanding antibiotic treatment criteria results in fewer antibiotic doses used over the course of an outbreak compared to current treatment practices (Figure 2). 


<div align="center">
  <img src="images/CholeraEpiResults.png" alt="Figure2" width="500"/>
</div>

Figure 2: Plot of population-level impact of expanded antibiotic treatment guidelines. Each plot shows the proportion of infections averted by the proportion of additional antibiotic doses used by proportion of all moderate infections who receive treatment (5%, 25%, 50%, 75%, 100%) and R_e scenario (low, intermediate, high). The outcomes are split into three regions by the impact of expanded eligibility criteria: between the dashed line along the x-axis (expanding criteria averts no infections) and the diagonal dashed line (each addition dose of antibiotics deployed averts one infections), points colored yellow, in which expanded eligibility only has individual-level benefits; between the diagonal dashed line and the vertical dashed line (no additional doses are used to avert infections), points colored green, in which expanded eligibility results in each dose preventing more than one additional infection; and the region left of the vertical dashed line, points colored purple, in which expanded eligibility results in fewer doses used over the course of the outbreak than compared to current antibiotic treatment guidelines.







