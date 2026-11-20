# AERI Infrared Radiance Analysis and Quality Control: Barrow, 1 June 2022

## Project Overview
The Atmospheric Emitted Radiance Interferometer (AERI) is a ground-based Fourier-transform spectrometer operated by the U.S. DOE Atmospheric Radiation Measurement (ARM) programme. It measures downwelling infrared spectral radiance, which carries information on temperature, water vapour, trace gases and clouds.
This project screens and explores one ARM AERI file from Barrow, Alaska, covering 1 June 2022 (3880 records). The aim is to assess data quality and to see how the spectra behave under different sky conditions before any quantitative retrieval is attempted.

## Scientific Background
Gases absorb and emit infrared radiation through molecular rotational and vibrational transitions, and a layer in local thermodynamic equilibrium emits at its own temperature. Blackbody emission follows Planck's law:
```math
B_{\tilde{\nu}}(T)=\frac{2hc^{2}\tilde{\nu}^{3}}{\exp\left(\dfrac{hc\tilde{\nu}}{k_{B}T}\right)-1}
```
Inverting it for a measured radiance $`L_{\tilde{\nu}}`$ gives the brightness temperature:
```math
T_{b}=\frac{hc\tilde{\nu}/k_{B}}{\ln\left(1+\dfrac{2hc^{2}\tilde{\nu}^{3}}{L_{\tilde{\nu}}}\right)}
```
Here $`\tilde{\nu}`$ is wavenumber (m⁻¹), and $`h`$, $`c`$ and $`k_{B}`$ are the Planck constant, speed of light and Boltzmann constant. $`T_{b}`$ is the temperature of the blackbody that would produce the measured radiance at that wavenumber. It is not the air temperature, and equals it only for an optically thick, isothermal emitter. In the radiative-transfer picture, opaque channels see emission from nearby air, while transparent channels see higher, colder layers or cloud. This is background only: no radiative-transfer model is run here.

## Methodology
- **Processing:** radiance (`mean_rad`), wavenumber, instrument flags and surface meteorology are read from NetCDF. Fill values (-9999) become NaN and times are converted to UTC.
- **Units:** wavenumber in cm⁻¹ is multiplied by 100 to give m⁻¹, and radiance in mW m⁻² sr⁻¹ (cm⁻¹)⁻¹ by 10⁻⁵ to give SI. A Planck round trip and the Wien peak (1.961 T cm⁻¹) check the implementation.
- **Record status:** each record is *missing* (`missingDataFlag`), *hatch not open* (`hatchOpen` ≠ 1), *invalid spectrum* or *usable*.
- **Anomaly detection:** for 420–700 cm⁻¹, the RMS departure from a 7-record rolling median (across neighbouring records in time) is computed. With $`\sigma = 1.4826\,\mathrm{MAD}`$, spectra above median + 8σ are *suspect* (kept), and spectra above a fixed 10 mW m⁻² sr⁻¹ (cm⁻¹)⁻¹ are *glitches* (rejected).
- **Consistency check:** $`T_{b}`$(660–670 cm⁻¹) minus surface air temperature is rejected beyond 5σ of its median, for hatch-open records.
- **Hourly assessment:** record counts by UTC hour are compared with the *median hourly record count*, a relative measure and not the instrument's absolute expected completeness.

## Spectral Analysis
Six bands are examined, using band-mean radiance converted to $`T_{b}`$ at the mean wavenumber:
- **450–550 cm⁻¹:** far-infrared water-vapour rotation band.
- **660–670 cm⁻¹:** CO₂ ν₂ band centre. Strongly absorbing, so $`T_{b}`$ should lie near the near-surface air temperature, which motivates the consistency check.
- **800–1000 cm⁻¹:** atmospheric window, the channel most sensitive to cloud and the basis of the sky classification.
- **1030–1050 cm⁻¹:** ozone band, compared with the neighbouring continuum at **1100–1120 cm⁻¹**.
- **1500–1600 cm⁻¹:** water-vapour ν₂ band.

## Results and Visualizations
- **Record status (Fig. 6d):** of 3880 records, 3709 (95.6%) are usable, 164 (4.2%) have the hatch not open, 6 (0.2%) are invalid and 1 is missing. The invalid spectra cluster near 23.6 UTC (Fig. 6a).
- **Hourly availability (Figs. 4–5):** relative to a median of 163 records per hour, most hours are close to 100%. Usable records drop to about half in hour 17 UTC, where about half the hour has the hatch not open, and to about 75% in hour 23. Time steps are bimodal, near 18 s and 60 s.
- **Consistency check (Fig. 6b):** $`T_{b}`$(660–670) minus air temperature mostly lies between about -5 and +3 K, so the CO₂ band tracks surface air temperature to within a few kelvin, as expected for an opaque band. A few low-tail outliers are flagged.
- **Negative radiances (Fig. 6c):** under 0.5% of values in every 200 cm⁻¹ band, concentrated in 1000–1200 and 1600–1800 cm⁻¹ and absent below 800 cm⁻¹.
- **Time evolution (Fig. 3):** window radiance is low and steady ($`T_{b}`$ roughly 190–197 K) from about 6.5 to 13 UTC. It then rises to about 270 K, close to air temperature, by about 15.5 UTC and stays there until about 19 UTC. It varies strongly before 6.5 and after 19.5 UTC, while CO₂-band radiance changes little.
- **Composites (Figs. 1–2):** 1265 clear-sky candidates and 731 cloud-dominated candidates, leaving 1713 intermediate usable spectra. The cloud-dominated composite is smooth, with $`T_{b}`$ of about 269–273 K from 400 to 1800 cm⁻¹. The clear-sky composite sits near 276 K in the CO₂ band, falls to about 180 K in the window, and rises to roughly 230–245 K near the ozone band.

## Limitations and Future Work
- **Sky classification is exploratory.** "Clear-sky candidates" have window $`T_{b}`$ < 200 K and "cloud-dominated candidates" have $`T_{air}-T_{b}`$ < 10 K. These fixed single-channel thresholds depend on humidity, temperature profile and season, can be confused by thin cloud, haze or fog, and are not validated. No cloud optical depth, height, atmospheric profiles or gas concentrations are retrieved.
- The anomaly metric is elevated in the variable-sky periods (roughly 1–7 and 20–24 UTC; Figs. 3, 6a), so "suspect" spectra may reflect real sky variability as well as instrument noise. In Fig. 2, composite $`T_{b}`$ above 1400 cm⁻¹ shows spikes well above surface air temperature, which I have not yet traced.
- Band $`T_{b}`$ from band-mean radiance at the mean wavenumber is an unquantified approximation. The thresholds are empirical, `qc_time` is reported but not used for rejection, and the analysis covers one day.
- **Future work:** validate against ceilometer, lidar or cloud radar, and set clear-sky thresholds from line-by-line calculations (e.g. LBLRTM) driven by radiosondes. Further steps are spectral cloud indicators, extension to a full month, and comparison with ARM's own AERI quality flags.

## Skills Demonstrated
- Scientific Python (NumPy, pandas, SciPy, Matplotlib) for NetCDF instrument data, fill values, flags and time series.
- Spectral analysis with the Planck function, unit handling and robust statistics.
- Critical interpretation of observational data, including its limits.

Reference: Knuteson et al. (2004), *J. Atmos. Oceanic Technol.*, 21, 1763–1789 (AERI, Parts I–II).
