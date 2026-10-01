# Revised data and composition audit

- Frozen molecules: 1,035
- Tier counts: {'0': 177, '1': 562, '2': 125, '3': 171}
- Source counts: {'PPDB': 510, 'ECOTOX': 441, 'BPDB': 84}
- Nested endpoint invariant: True

## Composition proportions by tier

|   Tier |   agrochemical_flag:fungicide |   agrochemical_flag:herbicide |   agrochemical_flag:insecticide |   agrochemical_flag:other_agrochemical |   source:BPDB |   source:ECOTOX |   source:PPDB |   toxicity_type:Contact |   toxicity_type:Oral |   toxicity_type:Other |
|-------:|------------------------------:|------------------------------:|--------------------------------:|---------------------------------------:|--------------:|----------------:|--------------:|------------------------:|---------------------:|----------------------:|
|  0.000 |                         0.220 |                         0.458 |                           0.090 |                                  0.169 |         0.040 |           0.407 |         0.554 |                   0.616 |                0.254 |                 0.130 |
|  1.000 |                         0.237 |                         0.436 |                           0.082 |                                  0.238 |         0.109 |           0.377 |         0.514 |                   0.623 |                0.270 |                 0.107 |
|  2.000 |                         0.168 |                         0.200 |                           0.336 |                                  0.232 |         0.008 |           0.480 |         0.512 |                   0.488 |                0.376 |                 0.136 |
|  3.000 |                         0.023 |                         0.018 |                           0.719 |                                  0.135 |         0.088 |           0.567 |         0.345 |                   0.614 |                0.275 |                 0.111 |

## Tier-composition association tests

| Variable           |      Chi2 |   df |           P |   CramerV |        BH_q |
|:-------------------|----------:|-----:|------------:|----------:|------------:|
| source             |  38.9826  |    6 | 7.21405e-07 | 0.13723   | 1.08211e-06 |
| toxicity_type      |   9.06851 |    6 | 0.169762    | 0.0661885 | 0.169762    |
| herbicide          | 123.745   |    3 | 1.20443e-26 | 0.345775  | 3.6133e-26  |
| fungicide          |  40.1874  |    3 | 9.72382e-09 | 0.197049  | 1.94476e-08 |
| insecticide        | 338.806   |    3 | 3.95701e-73 | 0.572144  | 2.3742e-72  |
| other_agrochemical |  10.7679  |    3 | 0.0130495   | 0.101999  | 0.0156594   |