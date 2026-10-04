# Frozen pH ecology sensitivity

Dataset version: `EQ-PH-SHARED-v1.0.0`. Analysis version: `ph-group-linkage-v1.0.0`.

The archived workbook contributed 709 admitted pH measurements linked to ecology through 563 site-campaign-position groups. Physical specimen identity is unverified; the linkage table preserves the reported Trip 4 material relationships.
Group-mean pH ranged from 7.290 to 9.980.

The pipeline fits group-linked alpha and composition models. Its geographic diagnostic first residualizes group composition against campaign-by-position indicators, with or without group pH, averages residuals within site, then fits linear and quadratic route terms. Each reported R² describes its own residual-response matrix; the difference is a descriptive diagnostic. Legacy partial_r2 fields retain this statistic for compatibility.

Rows that are pending, depleted, date-quarantined, or quality-control-quarantined are absent from every model. Availability is non-random, so results are bounded to this fixed cohort.

The ecology and data papers use this same frozen source. Future corrections or recovered measurements require a successor version and a new comparison.
