# Biology behind the statistics (descriptive context)

- Status: `descriptive_biology_context_complete`

## Landform

- Sites per landform: {'sand dune': 44, 'saline pan': 9, 'desert oasis': 3, 'aeolian lake': 1, 'dune slack': 1, 'gravel': 1, 'oilspill': 1}
- Saline-pan sites by transect third: {'east': 5, 'west': 3, 'central': 1}
- Dune vs saline pan composition, adjustment none: pseudo-F 4.56, p 0.0013
- Dune vs saline pan composition, adjustment route_linear_quadratic: pseudo-F 1.74, p 0.0362
- Genera differing (BH q<0.05): 16 unadjusted, 0 route-adjusted
- Higher in saline pans (route-adjusted): 
- Higher in dunes (route-adjusted): 
- Route model R2, all_60_sites: 0.402
- Route model R2, sand_dune_sites: 0.342
- Route genera supported: 124 (all sites) vs 71 (dune sites); 70 supported in both with the same sign
- Climate-diversity correlations on dune sites only: all nine supported = False

## Compartment genus family (600 tests)

- Supported: 295 of 600; by contrast and direction {'Deep-Surface:shallow_subsurface': 33, 'Deep-Surface:surface': 39, 'Rhizosphere-Deep:root_adjacent': 53, 'Rhizosphere-Deep:shallow_subsurface': 61, 'Rhizosphere-Surface:root_adjacent': 54, 'Rhizosphere-Surface:surface': 55}
- shallow_subsurface vs surface: higher in shallow_subsurface: Polycyclovorans (+1.79), Opitutus (+1.16), Thermosipho (+1.02), Nordella (+0.92), Acidiferrimicrobium (+0.89), MND1 (+0.89), I-8 (+0.85), Aureibacillus (+0.84), Sandaracinus (+0.84), Caldilinea (+0.82)
  higher in surface: Planococcus (-1.75), Nibribacter (-1.48), Kineococcus (-1.36), Deinococcus (-1.30), Saccharibacillus (-1.27), Kineosporia (-1.20), Rufibacter (-1.20), Kocuria (-1.13), Massilia (-1.12), Arcticibacter (-1.04)
- root_adjacent vs surface: higher in root_adjacent: TM7a (+2.07), Pelagibacterium (+1.86), Arsenicitalea (+1.67), Ohtaekwangia (+1.60), Neorhizobium (+1.51), DSSF69 (+1.47), Litchfieldia (+1.43), Metabacillus (+1.41), Devosia (+1.31), MM2 (+1.30)
  higher in surface: Bradyrhizobium (-2.27), Deinococcus (-2.01), Oxalophagus (-1.75), Blastococcus (-1.73), Vallicoccus (-1.64), Rubrobacter (-1.58), Limnobacter (-1.55), Cystobacter (-1.54), Symbiobacterium (-1.48), Geodermatophilus (-1.36)
- root_adjacent vs shallow_subsurface: higher in root_adjacent: Pelagibacterium (+2.37), TM7a (+2.01), Planococcus (+1.99), Brevundimonas (+1.95), Devosia (+1.95), Pseudomonas (+1.87), Neorhizobium (+1.76), Nibribacter (+1.75), Rhodocista (+1.63), Pantoea (+1.63)
  higher in shallow_subsurface: Bradyrhizobium (-2.18), MND1 (-2.03), Nitrospira (-2.01), Thermosipho (-1.99), Gaiella (-1.94), Oxalophagus (-1.79), Symbiobacterium (-1.76), Rubrobacter (-1.63), Caldilinea (-1.58), Acidiferrimicrobium (-1.53)

## XRF elemental axis

- Positive loadings: ['Ca', 'Mg', 'Na', 'S', 'Cl', 'Fe', 'Ti']; negative: ['Si']
- Axis vs route position: rho +0.73 (p 4.2e-11); genera tracking the axis: 118
- Positive (evaporite/carbonate side): Polygonibacillus (+0.81), Escherichia-Shigella (+0.74), Sediminibacillus (+0.74), Aquibacillus (+0.74), Halalkalibacter (+0.73), Gracilibacillus (+0.72), Enterococcus (+0.71), Aeromonas (+0.69), Pseudalkalibacillus (+0.69), Crossiella (+0.65)
- Negative (quartz side): Steroidobacter (-0.71), Caldilinea (-0.67), Phenylobacterium (-0.65), Paraconexibacter (-0.65), Streptomyces (-0.65), Roseisolibacter (-0.63), Gemmatimonas (-0.63), Pirellula (-0.58), Blastocatella (-0.58), Dongia (-0.58)

## pH

- Sites 60; genera tracking site pH: 117 (104 also route-supported)
- Higher at higher pH: Polygonibacillus (+0.81), Sediminibacillus (+0.80), Halalkalibacter (+0.78), Gracilibacillus (+0.75), Halomonas (+0.75), Aquibacillus (+0.71), Aeromonas (+0.71), Ammoniphilus (+0.69), Escherichia-Shigella (+0.69), Truepera (+0.68)
- Higher at lower pH: Pirellula (-0.68), Caenimonas (-0.68), Roseisolibacter (-0.65), Gemmatimonas (-0.65), Steroidobacter (-0.64), Gaiella (-0.64), Dongia (-0.62), Adhaeribacter (-0.61), Sphingomonas (-0.60), Streptomyces (-0.60)

## Core genera per compartment (>= 90% of sites after subsampling to 12865 reads)

- surface: 82 core genera of 912 detected; mean occupancy 0.281; median 274 genera per site
- shallow_subsurface: 69 core genera of 937 detected; mean occupancy 0.282; median 284 genera per site
- root_adjacent: 85 core genera of 1026 detected; mean occupancy 0.283; median 299 genera per site
- Core shared by all three compartments: 59

## Permitted wording

- Landform contrasts are site-level and descriptive; saline pans sit at both ends of the route, so the route-adjusted contrast is the informative one.
- Compartment genus contrasts are paired within site and campaign; they describe which genera carry the compartment difference and do not identify a mechanism.
- The XRF axis is an elemental axis (Ca, Mg, Na, S, Cl, Fe, Ti positive; Si negative); do not call it salinity.
- Core counts describe occupancy after subsampling and do not test a hypothesis.
