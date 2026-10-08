Very small colonies may have a\* dominated by background pixels. Evidence in the data: in Pb, wells below about 1,850 px have a\* of 6-7 with SD 2-2.7 at every dose (a hinge fit of a\* on ln area finds a flat segment below 1,846 px, with 30% of Pb wells below it); in Cr, wells below about 1,100-1,800 px have a\* of 4-6 with SD under 3.6. Cu, Fe and Zn colonies are rarely that small (5th percentile of area: Cu 6,759 px, Fe 13,829 px).

To test how much this matters, well-images whose largest object is smaller than a minimum area (1,000, 2,000, 3,000 or 5,000 px) were dropped, and the dose models (strain, run and plate random effects, at fixed size) were refitted. 90%, 78%, 70% and 59% of well-images remain. The 2,000 px value is the main sensitivity setting; the others show how much the result depends on that choice. Zinc is excluded.

{{IMG:s8_minarea|Figure 18. Dose effect against area threshold (top), per-dose effect at fixed size by threshold (middle), and wells retained by dose (bottom).}}

{{T_S8_MODEL}}

*Table 18. Dose effect on a\* over the full dose range, at fixed size, by minimum colony area (0 = no filter).*

{{T_S8_CR}}

*Table 19. Chromium: a\* difference from 0 dose at fixed size, by threshold, and wells retained.*

{{T_S8_PB}}

*Table 20. Lead: a\* difference from 0 dose at fixed size, by threshold, and wells retained.*

{{T_S8_TOP}}

*Table 21. Strain-mean a\* at the top dose by threshold.*

**Reading.**
- **Cu and Fe are unaffected** (Cu -13.4 to -13.3; Fe +3.8 to +4.4). Their colonies are large at every dose.
- **Cr: the top-dose drop is not a small-colony artifact.** The effect at dose 1.2 is -9.1 with no filter and -10.7 with colonies of at least 5,000 px (98 wells, 74 strains remain). The effects at doses 0.2-0.8 do not change. The overall dose effect stays between -8.4 and -9.3.
- **Pb: the rise at doses 5-15 holds or strengthens** with a threshold (+4.9, +8.0 and +9.5 at doses 5, 10 and 15 at 2,000 px, against +4.4, +7.5 and +7.6 with no filter). At doses 20-30 only 242-265 wells remain at 2,000 px, 90 at 3,000 px and 8-17 at 5,000 px (572-657 at 1,000 px), so the inhibitory range of Pb cannot be tested. In those few colonies the fixed-size difference from control is between -2.6 and +1.2. The Pb overall dose effect is unstable across thresholds (-6.1, -3.7, -4.3 [singular fit], -6.0, +0.2) because it depends on how many top-dose wells are left, and should not be quoted.
- **A floor remains in Pb.** Among Pb colonies of at least 3,000 px at dose 30 (81 strains), strain-mean a\* is still 5.7 with SD 0.56 across strains. Area thresholds alone do not tell us whether this is background or a real loss of pigment.
- Fixed-size effects extrapolate a size-a\* slope fitted mostly on larger colonies; thresholds keep the comparison within the range where colonies of similar size exist, which is why the Pb inhibitory doses lose almost all their wells.
