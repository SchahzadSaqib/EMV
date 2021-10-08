## Supplementary tables

## Relationship of birth mode of prior delivery with single bacterial taxa 
## in women with one prior delivery compared to nulliparas.

GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Previous_delivery_method',
              compare.to = 'Nulliparous',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              select.by = 'Parity_0or1',
              select = 'Yes', 
              pdf = T,
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)


## Supplementary Table 6. Relationship of parity (continuous) 
## with single bacterial taxa.
CovariateTest_mod(species.table = EMMI_phy324,
                   meta = EMMI_phy324,
                   covariate = 'Parity',
                   readcount.cutoff = 0, 
                   outlier.cutoff = 3, 
                   p.cutoff = 0.05, 
                   min.prevalence = 0.05, 
                   min.abundance = 0.05, 
                   keep.result = F,
                   nonzero = F, 
                   relative = T, 
                   pdf = F)



## Supplementary Table 5. Relationship of time (months) since previous delivery 
## and current pregnancy at the time of sample retrieval 
## and vaginal microbiota (n=111).

CovariateTest_mod1(species.table = EMMI_phy324,
                   meta = EMMI_phy324,
                   covariate = 'Time_between_sample_and_prev_preg_mo',
                   readcount.cutoff = 0, 
                   select.by = 'Prev_preg_birth', 
                   select = 'yes',
                   outlier.cutoff = 3, 
                   p.cutoff = 0.05, 
                   min.prevalence = 0.05, 
                   min.abundance = 0.05, 
                   keep.result = F,
                   nonzero = F, 
                   relative = T, 
                   pdf = F)


# Supplementary Table 6. Relationship of gestational age 
# with single bacterial taxa.
# All samples (n=324)
CovariateTest_mod(species.table = EMMI_phy324,
                  meta = EMMI_phy324,
                  covariate = 'Ga_at_sample_weeks',
                  readcount.cutoff = 0, 
                  outlier.cutoff = 3, 
                  p.cutoff = 0.05, 
                  min.prevalence = 0.05, 
                  min.abundance = 0.05, 
                  keep.result = F,
                  nonzero = F, 
                  relative = T)


## Stratified analysis (Nulliparous with no pregnancies, 
## nulliparous with previous pregnancies, parous)
CovariateTest_mod(species.table = EMMI_phy324,
                  meta = EMMI_phy324,
                  covariate = 'Ga_at_sample_weeks',
                  readcount.cutoff = 0, 
                  group = 'GP3',
                  outlier.cutoff = 3, 
                  p.cutoff = 0.05, 
                  min.prevalence = 0.05, 
                  min.abundance = 0.05, 
                  keep.result = F,
                  nonzero = F, 
                  relative = T, 
                  pdf = F)


## Supplementary Table 7. Relationship of demographic variables 
## with single bacterial taxa.
## Level of education
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324,
              group = 'Higher_education',
              compare.to = 'No',
              readcount.cutoff = 0,
              outlier.cutoff = 3, 
              p.cutoff = 0.05,
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)


## Current smoking
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Smoking_now',
              compare.to = 'No',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)


## History of smoking
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Smoking_ever',
              compare.to = 'No',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)


## Sexpartners during lifetime
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Sexpartners_over3',
              compare.to = 'No',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)

#Infertility treatments in current pregnancy
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Infertility_treatment',
              compare.to = 'No',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)


## Infertility treatments ever
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Infert_treatments_ever', 
              compare.to = 'No',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)


## Intercourse <48 hours prior to sampling
GroupTest_mod(species.table = EMMI_phy324, 
              meta = EMMI_phy324, 
              group = 'Recent_intercourse',
              compare.to = 'No',
              readcount.cutoff = 0, 
              outlier.cutoff = 3, 
              p.cutoff = 0.05, 
              min.prevalence = 0.05, 
              min.abundance = 0.05, 
              label.direction = 1, 
              keep.result = F,
              nonzero = F, 
              relative = T, 
              show_quartz = F)