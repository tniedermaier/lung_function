# lung_function
Codes for estimations of lung function decline in the COSYCONET cohort of COPD patients.

We aimed at estimating the average annual decline in lung function measured by the forced expiratory volume in 1 second (FEV1). To do so, we used random slopes models. While the basic idea is simple, a number of options exist that can influence the results significantly and where no option is clearly preferred. Hereafter are some explanations of the analysis and the rationale behind the decisions made.

- The data set comprises 2741 participants with physician-confirmed COPD diagnosis.
- A wide range of diagnostic and interview data are available.

Lufu_Decline_cosy.txt:
- Lung function (FEV1) was measured at each visit, along with EQ-5D-3L and many other variables.
- GOLD grade was assigned according to FEV1. Tiffeneau index was considered in additional analyses but not primarily. Rationale: There was considerable "diagnostic instability", i.e., patients who would be classifed as COPD patients at one visit and as "pre-COPD" or PRISm or "healthy" in the subsequent visit (18 months later), or vice versa. This is a consequence of the cutoff values for Tiffeneau index (FEV1 %predicted / FVC %predicted <0.70) in combination with the variability in measurements (e.g., submaximal effort at spirometry, unrecognized mild exacerbation etc.).
- Exacerbations were coded according to the most severe exacerbation in the past 12 months, with 0=no exacerbation, 1=only mild exacerbations (could be handled by the patient himself by increasing the medication dose), 2=moderate exacerbations (required corticosteroids and/or a physician visit), 3=severe exacerbation (required hospitalization).
- Follow-up status including dropout reasons (if dropout occured) were recorded for all participants and used to characterize the cohort, e.g., average FEV1 decline stratified by dropout reason.
- To address regression to the mean, we used FEV1 measurements for decline estimations that were not used for the GOLD classification and vice versa (except for sensitivity analysis).

copd_progression_from_long_baseR.R:
- FEV1 trajectories were used to estimate the average annual (and 3-month) progression rates from one GOLD grade to the next more advanced GOLD grade.
