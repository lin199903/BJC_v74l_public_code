# Post hoc sensitivity addendum

This addendum was written after the first locked analysis run showed a weak positive
bulk association for the ACP5/MMP9 core. It is explicitly exploratory and does not
alter the locked primary analysis.

The sensitivity analysis asks whether the bulk association persists after adjustment
for general myeloid abundance. A leave-target-and-program-out pan-myeloid score is
defined as the mean within-cohort z score of `LST1`, `AIF1`, `FCER1G`, `TYROBP`,
`CTSS`, `LILRB1`, `CSF1R` and `SPI1`. All eight genes were chosen before inspecting
the adjusted estimates. Within each bulk cohort, ranks of SIGLEC15 and the ACP5/MMP9
score are separately residualized on the rank of this pan-myeloid score; their
residual correlation is reported with a 10,000-resample paired bootstrap CI. The
three partial correlations are combined using the same REML--Hartung--Knapp method
as the locked analysis.

This analysis can diagnose composition dependence but cannot identify a causal
cell state or substitute for measured tumour purity or same-sample protein assays.
