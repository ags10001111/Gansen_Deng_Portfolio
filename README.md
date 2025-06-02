# 👋 Hi, I'm Gansen!

I am a **PhD in Statistics**. My PhD research focused on developing statistical learning methods for self-reported data, with applications to chronic pain research. My academic journey has equipped me with a strong foundation in **statistical theory**, **machine learning**, and **data analysis**—skills I’m excited to apply to real-world challenges in industry.

Throughout my research, I’ve worked on:
- 🧠 **Clustering** and **metric learning**
- 🏷️ **Multi-label classification**
- 🔬 Simulation studies and modeling with complex error structures
- 📊 Working with noisy, self-reported, and high-dimensional data

These experiences have sharpened my analytical thinking and programming proficiency in languages such as **R** and **Python**.

I’m passionate about turning data into actionable knowledge and building models that support better decision-making.

---

📁 This repository is a portfolio of my:
- Research projects
- Data science explorations
- Learning progress and tools
- Technical skillset

📄 [View My CV (PDF)](CV_GansenDeng.pdf)

# Table of Contents


# Published Papers

## A real-time and interactive web-based platform for visualizing and analyzing COVID-19 in Canada (2020)

- 🌐 [Website Link](https://covid-19-canada.uwo.ca/)
- 📄 [Paper Link](https://www.ccsenet.org/journal/index.php/ijsp/article/view/0/43346)

This website was developed during the COVID-19 pandemic to monitor real-time changes in the spread of the virus across Canada. It was created by our research group, the GW Data Science Research Group (GW-DSRG). My primary responsibility was to design and implement the real-time data visualizations for the provinces of [Ontario](https://covid-19-canada.uwo.ca/en/ontario.html) ([codes](COVID19_Website/IR_ON.R)) and [Alberta](https://covid-19-canada.uwo.ca/en/alberta.html) ([codes](COVID19_Website/IR_AL.R)). The interactive plots were built using the R package **plotly**.

## Epilepsy-associated death in Southwestern Ontario (2023)
- 📄 [Paper Link](https://onlinelibrary.wiley.com/doi/full/10.1111/bpa.13121)

**Comparison analysis** ([codes](Epilepsy_Study/Comparison.R)) was performed on this dataset. Continuous variables were compared between groups using the **Mann–Whitney U test**, while discrete variables were analyzed using **Fisher's exact test**. **Logistic regression** ([codes](Epilepsy_Study/GLM.R)) was also used to model the cause of death, with **LASSO regularization** applied for variable selection.


## Longitudinal analysis of mucosa-associated invariant T cells in sepsis study (2023)
- 📄 [Paper Link](https://pubmed.ncbi.nlm.nih.gov/36604951/)
- 🖼️ [Poster](Sepsis_Study/Sepsis_Poster.pdf)

![Sepsis Poster](Sepsis_Study/Sepsis_Poster.jpg)

Data on MAIT cells was collected via flow cytometry from blood samples taken across six time points from three cohorts: *septic, non-septic*, and *healthy*. Three types of analyses were conducted:

1. **Comparative analysis** ([codes](Sepsis_Study/Comparison.R)): T-tests and Mann–Whitney U tests were used to identify variables that differed significantly between cohorts.  
2. **Longitudinal analysis** ([codes](Sepsis_Study/Longitudinal_analysis.R)): Generalized estimating equations (GEE) were applied to assess significant trends in the variables over time.  
3. **Correlation analysis** ([codes](Sepsis_Study/Correlation.R)): Logistic regression was used to explore associations between variables, with variable selection performed using the likelihood ratio test.

# PhD Thesis 

- 📄 [Thesis Link](https://ir.lib.uwo.ca/etd/10805/) (Statistical Learning Methods for Challenges arised from Self-Reported Data)
- Goal: Develop data-driven methods to assist clinicians in managing chronic pain

## Phenotyping Chronic Pain Patients using a Mechanism-Based Framework
- **Objective:**
  - Cluster patients into phenotypes using the routinely clinically collected variables
  - Check the ability of these variables to differentiate the phenotypes & mechanism-based classification (*non-nociplastic, nociplastic, mixed*)

- **Methods:**
  - The **latent class analysis (LCA)** ([codes](Phenotyping_CP_Patients/LCA.R)) is applied to cluster CP data with 13 variables and 198 patients
  - The optimal number of clusters in LCA is determined using **Bayes' information criteria (BIC)**
  - The **chi-square independence test** is used to check whether the features are significantly different across clusters
  - The feature importance is checked using **random forest** models by treating the clustering labels/mechanism-based classification as the responses. The **permutation importance** is used to quantify the importance of each feature

## A Novel Distance Metric for Clustering Questionnaire Data
- **Objective:**
  - Develop a distance metric tailored for mixed-type data comprising both continuous and categorical variables
  - Design a metric that incorporates the unique characteristics (the rating criteria vary between individuals but remain consistent within each individual) of questionnaire-based data to mitigate subjective bias and more accurately capture similarity between subjects

- **Notation:**

<p align="center">
  <img src="Novel_Distance/Notation.png" alt="Notation Table">
</p>

- **The Proposed Distance:** ([codes](Novel_Distance/Distance_Definition.R))

The distance between two subjects $x_i$ and $x_j$ is proposed as:

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\dpi{150}&space;\color{White}d(\mathbf{x}_i,&space;\mathbf{x}_j)&space;=&space;\frac{1}{p&space;&plus;&space;s}&space;\left(&space;\sum_{k=1}^{m_1&space;&plus;&space;m_2}&space;\frac{|x_{ik}&space;-&space;x_{jk}|}{R_k}&space;&plus;&space;\sum_{t=1}^{s}&space;w_t&space;\sin&space;\left(\frac{\arccos(r_{ij,t})}{2}&space;\right)&space;&plus;&space;\sum_{k&space;=&space;m_1&space;&plus;&space;m_2&space;&plus;&space;1}^{p}&space;\delta_{\text{cat}}(x_{ik},&space;x_{jk})&space;\right)" alt="Equation">
</p>

where:
- $\delta_{\text{cat}}(x_{ik}, x_{jk})$ is a co-occurrence-based distance [Ahmad and Dey, 2007](https://www.sciencedirect.com/science/article/abs/pii/S0169023X0700050X) for categorical variables.
- $w_t$ is the standard deviation of the lower triangular elements of the correlation distance matrix for the $t$-th group of self-reported (SR) variables, where the correlation distance is defined as $\sin\left( \frac{\arccos(r_{ij,t})}{2} \right)$

---

- **Selected Simulation Results:** ([codes](Novel_Distance/Simulation_Study.R))

The proposed distance achieves the highest ARI compared to all other distance metrics across all scenarios.

<p align="center">
  <img src="Novel_Distance/ARI_Normal.png" alt="ARI_Normal" width="800">
</p>

## Chronic Pain Patient Clustering by Accommodating Self-report Questionnaires and Interpretability
- **Objective:**
  - Cluster the CP patients based on the proposed questionnaire distance to capture the unique characteristics of the self-reported (SR) variables
  - Make clustering results more interpretable

- **Methods:**
  - Hierarchical clustering (HC) is employed with complete linkage, and the optimal number of clusters is determined by maximizing the Silhouette score
  - Interpretable Clustering via Optimal Trees (ICOT) [Bertsimas et al., 2021](https://link.springer.com/article/10.1007/S10994-020-05896-2) is applied to generate an interpretable clustering result with a tree structure ([codes](Interpretable_Clustering/icot.py))
  - The cluster center is computed for each cluster to represent its characteristics
  - Feature importance from the random forest model is used to evaluate the importance of each feature

- **ICOT Results:**

<p align="center">
  <img src="Interpretable_Clustering/ICOT_Result.png" alt="ARI_Normal" width="800">
</p>

## Semi-supervised Clustering of Self-reported Data using Active Learning

## Clinical Characteristics of Myofascial Trigger Points
- **Objective:**
  - Identify key MTrP characteristics crucial for clinical diagnosis and research  
  - Investigate the correlations among MTrPs, pain-pressure thresholds (PPT), and Michigan Body Map

- **Methods:**
  - The canonical correlation analysis (CCA) is used to investigate the relationships among the three sets of variables

## Clinical Prediction of Nociplastic Pain using Combined Criteria

- **Objective:**
  - Assess the predictive ability of collected variables for diagnosing nociplastic pain  
  - Identify optimal threshold values for those variables to establish simple diagnostic rules
 
- **Methods:**
  - Three prediction models are considered: logistic regression, random forest (RF), and support vector machine (SVM). The data is randomly split into 70% for training and 30% for testing. This train/test split is repeated 100 times, and the models are evaluated based on the mean prediction accuracy on the test data
  - Use the best prediction model to determine the optimal cutoff for each predictor that yields the highest mean prediction accuracy

# Other Projects







