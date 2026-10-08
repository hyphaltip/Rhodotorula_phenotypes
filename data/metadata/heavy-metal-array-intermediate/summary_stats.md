# Summary statistics: heavy-metal-array-intermediate

Total rows: 517,371; columns in union: 179

| Metal | rows | cols in file | runs | strains (distinct strain_id) | NULL strain_id | conc levels | conc min-max | first date | last date |
|---|---|---|---|---|---|---|---|---|---|
| Chromium | 168,793 | 147 | 5 | 334 | 5,923 | 7 | 0-1.2 | 2026-02-01 | 2026-02-06 |
| Copper | 155,470 | 147 | 5 | 334 | 3,489 | 7 | 0-30 | 2026-02-20 | 2026-02-25 |
| Iron | 66,506 | 147 | 3 | 243 | 0 | 7 | 0-30 | 2026-05-16 | 2026-05-20 |
| Lead | 113,386 | 179 | 5 | 333 | 1,509 | 7 | 0-30 | 2026-05-09 | 2026-05-14 |
| Zinc | 13,216 | 178 | 2 | 170 | 0 | 7 | 0-30 | 2026-05-01 | 2026-05-06 |

## Runs per metal

- Chromium d000320: 47,207
- Chromium d000321: 42,083
- Chromium d000322: 42,230
- Chromium d000323: 32,431
- Chromium d000324: 4,842
- Copper d000353: 43,415
- Copper d000354: 38,556
- Copper d000355: 37,990
- Copper d000356: 32,526
- Copper d000357: 2,983
- Iron d000406: 32,877
- Iron d000407: 27,618
- Iron d000408: 6,011
- Lead d000399: 34,733
- Lead d000400: 26,963
- Lead d000401: 26,099
- Lead d000402: 24,761
- Lead d000403: 830
- Zinc d000388: 6,458
- Zinc d000390: 6,758

## Missing values in key columns (all metals)

| column | NULL count |
|---|---|
| Shape_Area | 0 |
| ColorLab_L*GeoMedian | 0 |
| Intensity_MeanIntensity | 0 |
| strain_id | 10,921 |
| Concentration | 0 |
| capture_datetime | 0 |

## Columns only in Pb and Zn (32)

`GridLinReg_ColB`, `GridLinReg_ColM`, `GridLinReg_PredCC`, `GridLinReg_PredRR`, `GridLinReg_ResidualError`, `GridLinReg_RowB`, `GridLinReg_RowM`, `OrientZones_OutwardRotationConsistency-Mask-Dense`, `OrientZones_OutwardRotationConsistency-Mask-Overall`, `OrientZones_OutwardRotationConsistency-Mask-Sparse`, `OrientZones_OutwardRotationNet-Mask-Dense`, `OrientZones_OutwardRotationNet-Mask-Overall`, `OrientZones_OutwardRotationNet-Mask-Sparse`, `OrientZones_OutwardRotationRate-Mask-Dense`, `OrientZones_OutwardRotationRate-Mask-Overall`, `OrientZones_OutwardRotationRate-Mask-Sparse`, `OrientZones_OutwardRotationSustainedPeak-Mask-Dense`, `OrientZones_OutwardRotationSustainedPeak-Mask-Overall`, `OrientZones_OutwardRotationSustainedPeak-Mask-Sparse`, `Size_Area`, `Size_IntegratedIntensity`, `SymZones_CoreArea`, `SymZones_CoreEndRadius`, `SymZones_CoreRadius`, `SymZones_DenseArea`, `SymZones_DenseEndRadius`, `SymZones_MaxExpansion`, `SymZones_MeanExpansion`, `SymZones_SparseArea`, `SymZones_SparseEndRadius`, `SymZones_SymmetricRadius`, `index`
