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

## Epilepsy-associated death in the Southwestern Ontario (2023)
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

## Phenotyping Chronic Pain Patients using a Mechanism-Based Framework

## A Novel Distance Metric for Clustering Questionnaire Data

## Chronic Pain Patient Clustering by Accommodating Self-report Questionnaires and Interpretability

## Semi-supervised Clustering of Self-reported Data using Active Learning

## Clinical Characteristics of Myofascial Trigger Points

## Clinical Prediction of Nociplastic Pain using Combined Criteria



# Other Projects







