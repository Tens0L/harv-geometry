# HARV nose geometry

Nose profiles of the High-Speed Army Reference Vehicle (HARV), following
ARL-TN-1128 (Vasile et al., DEVCOM ARL, 2022), Appendix A.
<https://apps.dtic.mil/sti/pdfs/AD1176136.pdf>

All values are in calibers (base diameter D = 1).

| Parameter | Value |
|---|---|
| Nose length L | 5 (fineness ratio 5) |
| Base radius r_n | 0.5 |
| Spherical tip radius r_b | 0.05 |

| Profile | x_t | r_t | Note |
|---|---|---|---|
| Von Kármán (Sears–Haack, C=0) | 0.039103 | 0.048798 | tangent-blunted |
| Cone (power series n=1) | 0.045475 | 0.049795 | half-angle 5.192° |

## Files
- `harv_nose_vonkarman.csv`, `harv_nose_cone.csv` — (x, r) profile points
- `harv_nose.py` — generator (`python3 harv_nose.py --D 0.1` to scale)
- `harv_nose.png` — profile comparison plot
- `SUMMARY.md` — Japanese summary of ARL-TN-1128
- `harv_nose_excel.xlsx` — Excel version built only with worksheet formulas
- `EXCEL_GUIDE.md` — how to build the profiles in Excel (Japanese)
