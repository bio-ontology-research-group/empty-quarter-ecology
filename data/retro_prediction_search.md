# Retro-prediction search — summary and negative finding

**Task:** find a published 16S experiment with (a) a salinity manipulation or (b)
a desert-soil rewetting, where the digital twin's sign prediction can be
retro-tested.  Conclusion: **no clean match exists in the public literature
within the twin's validated regime.**  This memo documents the search so that
the finding can be cited honestly if the question is raised during review.

## Twin's validated regime (to match a candidate experiment against)

The twin identifies credibly-signed Shannon responses for:

1. **Acute rainfall pulse** at 5--30 mm natural intensity, 7--14 d lag,
   surface compartment only (Hill $K \approx 9.5$ mm, $B_{\max} \approx -0.54\sigma$).
2. **Salinity reduction** at the 1--2 SD scale of the Rub' al-Khali Trip-5
   XRF gradient ($\Delta$S roughly 0.05--0.25\% dry mass).
3. **Combined reclamation** do($S-2$ SD, $P+1$ SD): all compartments
   credibly positive, $+0.38$ to $+0.55\sigma$ Shannon.

A useful retro-prediction target has to land inside one of these regimes.

## Candidates examined

| Paper | Regime | V-region | Per-sample Shannon | Public data | Verdict |
|------|--------|----------|--------------------|-------------|---------|
| Atacama simulated rainfall (PMC10537920; DRR465014--087) | 50 mm eq., 30 d microcosm | V4 | increased over time | yes (DDBJ) | **out of regime** — 50 mm is far above the Hill saturation; Atacama rock pavement, not sand |
| Namib wetting frequency/intensity (Frossard 2015, srep12263) | fog/light/heavy pulses, 36 d | — | T-RFLP only | — | **wrong assay** (no Shannon from OTU table) |
| Junggar precipitation gradient (PMC6991129; SRR10020091--105) | $\pm 30\%$, $\pm 60\%$ over 2.5 y | V3--V4 | chronic, single endpoint | yes | **wrong regime** — chronic not acute-pulse; reports positive Shannon--rainfall relationship, opposite of twin's acute prediction but not contradictory (different mechanism) |
| Thar desert surveys (PMC8945486, PMC4729749) | observational | V3--V4 or legacy | aggregated only | limited | too small / not a manipulation |
| Saline-alkali subsurface pipe desalination | field | V3--V4 | per-treatment only | partial | **wrong biome** (agricultural, not desert sand) |
| Date palm root saline-irrigation (nature.com/41598-022-16869-x) | saline irrigation | — | — | — | out of biome |

## Why the regime is so narrow

The twin's diagnostic signature is a specifically-sand-desert, specifically-acute,
specifically-7--14-day-lag response in the *surface* compartment, driven by a
large pool of *conditionally rare genera* (408 genera, prevalence $<20\%$,
max $\geq 1\%$) whose bloom flattens Shannon after natural rainfall.  Every
published hyperarid rewetting experiment in the literature uses a
\emph{microcosm} design with $\geq 50$~mm water equivalents and reports
\emph{Shannon increases} at steady state --- consistent with reactivating a
dormant pool from near-zero diversity, but measuring a different quantity.
The $5$--$30$~mm natural-pulse regime that identifies our signal has not
been reproduced experimentally in any desert.

## What *is* available as weak sign-retro-prediction

The existing Atacama EC--Shannon correlation ($\rho = -0.38$, $p = 0.012$;
Fig.~3h) already constitutes a sign-level retro-prediction: the twin's
mediation identifies salinity as a negative driver of Shannon, and the
Atacama signal matches direction and magnitude on an independent continent.
This is documented in the cross-desert Results section but not currently
framed as "retro-prediction" --- we could reframe the sentence.

## Recommendation

1. **Do not add a new retro-prediction section.**  No candidate strong
   enough to carry the claim.
2. Instead, **reframe the existing Atacama EC--Shannon replication** in the
   cross-desert section explicitly as "the twin's sign prediction for the
   salinity$\to$Shannon edge is replicated on an independent continent",
   costing one sentence.
3. If a reviewer asks for experimental validation, acknowledge the gap:
   the twin's validated regime does not overlap any published hyperarid-sand
   manipulation, and purpose-designed field trials are proposed as the next
   step.

## Search log (for reference)

- Google Scholar / WebSearch queries, 2026-04-19, night run.
- Keywords: "desert soil 16S amplicon Shannon diversity salinity electrical
  conductivity gradient", "desert arid soil rewetting experiment 16S rRNA
  bacterial diversity Shannon decline pulse precipitation simulation",
  "saline soil reclamation desalination 16S amplicon bacterial diversity
  Shannon desert arid field experiment before after".
- 30 candidates scanned, none inside the twin's validated regime.
