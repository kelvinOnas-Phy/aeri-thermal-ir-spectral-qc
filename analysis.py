import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import netcdf_file
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
h = 6.62607015e-34
c = 2.99792458e8
k = 1.380649e-23
def brightness_temperature(radiance_mw, wavenumber_cm):
nu = wavenumber_cm * 100.0
B = np.where(radiance_mw > 0, radiance_mw * 1e-5, np.nan)
return (h * c * nu / k) / np.log(1.0 + 2.0 * h * c**2 * nu**3 / B)
# ===============================================================
# Section 1: data inspection and preparation
# ===============================================================
nc = netcdf_file("data/ARM_file.nc", mmap=False)
wnum = np.array(nc.variables["wnum"][:], dtype=float)
radiance = np.array(nc.variables["mean_rad"][:], dtype=float)
radiance[radiance == -9999] = np.nan
base_time = int(nc.variables["base_time"][()])
time_offset = np.array(nc.variables["time_offset"][:], dtype=float)
time_utc = pd.to_datetime(base_time + time_offset, unit="s")
obs = pd.DataFrame(index=time_utc)
for name, column in [("missingDataFlag", "missing"), ("hatchOpen", "hatch"),
("qc_time", "qc_time"), ("outsideAirTemp", "air_temp"),
("atmosphericPressure", "pressure"),
("atmosphericRelativeHumidity", "rel_humidity")]:
obs[column] = np.array(nc.variables[name][:], dtype=float)
obs = obs.replace(-9999, np.nan)
print("Site:", nc._attributes["facility_id"].decode(), "| datastream:",
nc._attributes["datastream"].decode())
print("First spectrum:", obs.index[0], "| last spectrum:", obs.index[-1])
print("Spectra:", radiance.shape[0], "| spectral points:", radiance.shape[1])
print("Wavenumber range: %.2f to %.2f cm-1, spacing %.3f cm-1" % (wnum.min(),
wnum.max(), np.diff(wnum).mean()))
print("NaN or fill values in radiance:", int(np.isnan(radiance).sum()))
print("Radiance range: %.1f to %.1f mW/(m2 sr cm-1)" % (np.nanmin(radiance),
np.nanmax(radiance)))
# ===============================================================
# Section 2: quality control, hourly
# ===============================================================
band_co2 = (wnum >= 660) & (wnum <= 670)
band_window = (wnum >= 800) & (wnum <= 1000)
band_signal = (wnum >= 420) & (wnum <= 700)
obs["bt_co2"] = brightness_temperature(radiance[:, band_co2].mean(axis=1),
wnum[band_co2].mean())
obs["bt_window"] = brightness_temperature(radiance[:, band_window].mean(axis=1),
wnum[band_window].mean())
obs["n_negative"] = (radiance < 0).sum(axis=1)
obs["flag_missing"] = obs["missing"] == 1
obs["flag_hatch"] = obs["hatch"] != 1
obs["flag_time"] = obs["qc_time"] != 0
hatch_open = obs["hatch"] == 1
signal = pd.DataFrame(radiance[:, band_signal])
local_median = signal.rolling(7, center=True, min_periods=3).median()
obs["anomaly"] = np.sqrt(((signal - local_median) ** 2).mean(axis=1)).values
a = obs.loc[hatch_open, "anomaly"]
anomaly_median = a.median()
anomaly_sigma = 1.4826 * (a - anomaly_median).abs().median()
suspect_limit = anomaly_median + 8 * anomaly_sigma
glitch_limit = 10.0
obs["flag_suspect"] = obs["anomaly"] > suspect_limit
obs["flag_glitch"] = obs["anomaly"] > glitch_limit
obs["bt_diff"] = obs["bt_co2"] - obs["air_temp"]
d = obs.loc[hatch_open, "bt_diff"]
diff_median = d.median()
diff_sigma = 1.4826 * (d - diff_median).abs().median()
obs["flag_consistency"] = hatch_open & ((obs["bt_diff"] - diff_median).abs() > 5
* diff_sigma)
obs["flag_invalid"] = obs["flag_glitch"] | obs["flag_consistency"]
obs["category"] = np.select(
[obs["flag_missing"], obs["flag_hatch"], obs["flag_invalid"]],
["missing", "hatch not open", "invalid spectrum"], default="usable")
obs["usable"] = obs["category"] == "usable"
hour = obs.index.hour
categories = ["usable", "invalid spectrum", "hatch not open", "missing"]
hourly = obs.groupby([hour,
"category"]).size().unstack(fill_value=0).reindex(columns=categories,
fill_value=0)
hourly["records"] = hourly.sum(axis=1)
reference = hourly["records"].median()
hourly["present_pct"] = 100 * hourly["records"] / reference
hourly["usable_pct"] = 100 * hourly["usable"] / reference
hourly["suspect"] = obs.groupby(hour)["flag_suspect"].sum()
usable_obs = obs[obs["usable"]]
hourly["neg_per_spectrum"] = usable_obs.groupby(usable_obs.index.hour)
["n_negative"].mean()
hourly.index.name = "hour_utc"
hourly.round(1).to_csv("results/hourly_qc.csv")
gap = obs.index.to_series().diff().dt.total_seconds().dropna()
neg_fraction_usable = (radiance[obs["usable"].values] < 0).mean()
print("\n--- Quality control ---")
print("Reference records per hour (median): %.0f" % reference)
print(hourly.round(1).to_string())
print("\nRecord categories:")
print(obs["category"].value_counts().to_string())
print("Hatch flag values:", obs["hatch"].value_counts().to_dict())
print("Records with qc_time flag:", int(obs["flag_time"].sum()), "at",
list(obs.index[obs["flag_time"]]))
print("Time step: median %.0f s, max %.0f s, steps over 120 s: %d" %
(gap.median(), gap.max(), (gap > 120).sum()))
print("Steps of 50-70 s: %d (spectra / 8 = %.0f)" % (((gap >= 50) & (gap <=
70)).sum(), len(obs) / 8))
print("Anomaly metric: median %.2f, sigma %.2f, suspect limit %.2f, glitch limit
%.0f" % (anomaly_median, anomaly_sigma, suspect_limit, glitch_limit))
print("Suspect spectra (kept):", int(obs["flag_suspect"].sum()), "| glitch:",
int(obs["flag_glitch"].sum()),
"| failing BT-CO2 vs air temperature check:",
int(obs["flag_consistency"].sum()))
print("BT(660-670) minus air temperature, hatch open: median %.2f K, sigma %.2f
K" % (diff_median, diff_sigma))
print("Invalid spectra:", list(obs.index[obs["flag_invalid"] &
hatch_open].strftime("%H:%M:%S")))
print("Negative radiance values in usable spectra: %.3f %%" % (100 *
neg_fraction_usable))
print("Hatch-closed spectra mean window radiance: %.1f vs hatch-open: %.1f" % (
np.nanmean(radiance[(obs["hatch"] == 0).values][:, band_window]),
np.nanmean(radiance[obs["usable"].values][:, band_window])))
# ===============================================================
# Section 3 and 4: radiance and brightness temperature
# ===============================================================
usable = obs["usable"].values
clear = usable & (obs["bt_window"] < 200).values
cloud_dominated = usable & ((obs["air_temp"] - obs["bt_window"]) < 10).values
rad_clear = radiance[clear].mean(axis=0)
rad_cloud = radiance[cloud_dominated].mean(axis=0)
bt_clear = brightness_temperature(rad_clear, wnum)
bt_cloud = brightness_temperature(rad_cloud, wnum)
print("\n--- Planck function checks ---")
T_test, nu_test = 300.0, 1000.0
B_test = (2 * h * c**2 * (nu_test * 100) ** 3 / (np.exp(h * c * nu_test * 100 /
(k * T_test)) - 1)) / 1e-5
print("Blackbody 300 K at 1000 cm-1: %.3f mW/(m2 sr cm-1); recovered temperature:
%.6f K" % (
B_test, brightness_temperature(B_test, nu_test)))
nu_grid = np.linspace(10, 3000, 100000)
B_grid = 2 * h * c**2 * (nu_grid * 100) ** 3 / (np.exp(h * c * nu_grid * 100 / (k
* 273.15)) - 1)
print("Peak of Planck curve at 273.15 K: %.1f cm-1 (Wien: 1.961 x T = %.1f cm-1)"
% (nu_grid[np.argmax(B_grid)], 1.961 * 273.15))
print("Overall usable: window radiance max %.1f, min of band means %.1f" % (
np.nanmax(radiance[usable][:, band_window].mean(axis=1)),
np.nanmin(radiance[usable][:, band_window].mean(axis=1))))
print("\n--- Sky classification from window brightness temperature ---")
print("Usable spectra:", usable.sum())
print("Clear-sky candidates (window BT < 200 K):", clear.sum())
print("Cloud-dominated candidates (window BT within 10 K of air temperature):",
cloud_dominated.sum())
print("Neither (intermediate):", usable.sum() - clear.sum() -
cloud_dominated.sum())
for limit in [195, 200, 205, 210, 220]:
print(" window BT < %d K: %d spectra" % (limit, (usable & (obs["bt_window"]
< limit).values).sum()))
classes = pd.DataFrame({"clear": clear, "cloud_dominated": cloud_dominated},
index=obs.index)
print(classes.groupby(classes.index.hour).sum().T.to_string())
# ===============================================================
# Section 5: clear-sky spectral features
# ===============================================================
bands = [("Far-IR water vapour rotation", 450, 550), ("CO2 band centre", 660,
670),
("Atmospheric window", 800, 1000), ("Ozone band", 1030, 1050),
("Continuum beside ozone", 1100, 1120), ("Water vapour nu2 band", 1500,
1600)]
print("\n--- Band brightness temperatures of the composite spectra ---")
print("%-30s %9s %9s %9s" % ("band (cm-1)", "clear K", "cloud K", "diff K"))
for label, lo, hi in bands:
m = (wnum >= lo) & (wnum <= hi)
bt_c = brightness_temperature(rad_clear[m].mean(), wnum[m].mean())
bt_k = brightness_temperature(rad_cloud[m].mean(), wnum[m].mean())
print("%-22s %4d-%4d %9.1f %9.1f %9.1f" % (label, lo, hi, bt_c, bt_k, bt_k -
bt_c))
ozone = (wnum >= 1030) & (wnum <= 1050)
beside = (wnum >= 1100) & (wnum <= 1120)
print("Clear-sky radiance, ozone band minus neighbouring continuum: %.2f mW/(m2
sr cm-1)" % (
rad_clear[ozone].mean() - rad_clear[beside].mean()))
print("Mean air temperature, clear / cloud-dominated: %.1f K / %.1f K" % (
obs.loc[clear, "air_temp"].mean(), obs.loc[cloud_dominated,
"air_temp"].mean()))
limit = obs.loc[clear, "air_temp"].max() + 5
print("Channels with composite BT above warmest clear-sky air temperature + 5 K
(%.1f K):" % limit)
for lo, hi in [(400, 1300), (1300, 1500), (1500, 1800)]:
m = (wnum >= lo) & (wnum < hi)
print(" %d-%d cm-1: clear %d of %d, cloud-dominated %d" % (lo, hi,
(bt_clear[m] > limit).sum(), m.sum(), (bt_cloud[m] > limit).sum()))
print("Correlation of BT(660-670) with air temperature (usable): %.2f" %
obs.loc[obs["usable"], ["bt_co2", "air_temp"]].corr().iloc[0, 1])
# ===============================================================
# Section 6: what the file says about clouds
# ===============================================================
cloud_names = [n for n in nc.variables if any(w in n.lower() for w in ["cloud",
"optical", "tau", "cbh", "lwp"])]
print("\n--- Cloud variables found in file ---")
print(cloud_names if cloud_names else "none")
# ===============================================================
# Section 7: figures
# ===============================================================
hours = obs.index.hour + obs.index.minute / 60 + obs.index.second / 3600
rad_label = r"Radiance [mW m$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$]"
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(wnum, rad_cloud, color="tab:orange", lw=0.8, label="Cloud-dominated
(n=%d)" % cloud_dominated.sum())
ax.plot(wnum, rad_clear, color="tab:blue", lw=0.8, label="Clear-sky candidates
(n=%d)" % clear.sum())
ax.set_xlabel(r"Wavenumber [cm$^{-1}$]")
ax.set_ylabel(rad_label)
ax.set_title("Mean downwelling radiance, Barrow, 1 June 2022")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("figures/fig1_radiance_spectra.png", dpi=200)
plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(wnum, bt_cloud, color="tab:orange", lw=0.8, label="Cloud-dominated")
ax.plot(wnum, bt_clear, color="tab:blue", lw=0.8, label="Clear-sky candidates")
ax.axhline(obs.loc[clear, "air_temp"].mean(), color="grey", ls="--", lw=0.8,
label="Mean surface air temperature")
ax.axvline(667, color="k", ls=":", lw=0.8)
ax.axvline(1042, color="k", ls=":", lw=0.8)
ax.text(672, 178, r"CO$_2$ 667", fontsize=8)
ax.text(1047, 178, r"O$_3$ 1042", fontsize=8)
ax.set_ylim(170, 300)
ax.set_xlabel(r"Wavenumber [cm$^{-1}$]")
ax.set_ylabel("Brightness temperature [K]")
ax.set_title("Brightness temperature of the composite spectra")
ax.legend(loc="lower right", fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("figures/fig2_brightness_temperature.png", dpi=200)
plt.close(fig)
rad_5min = pd.DataFrame(radiance[usable],
index=obs.index[usable]).resample("5min").mean()
full_day = pd.date_range(obs.index[0].normalize(), periods=288, freq="5min")
rad_5min = rad_5min.reindex(full_day)
fig, axes = plt.subplots(3, 2, figsize=(9, 10), gridspec_kw={"width_ratios": [40,
1]})
ax = axes[:, 0]
ax[1].sharex(ax[0])
ax[2].sharex(ax[0])
axes[0, 1].axis("off")
axes[1, 1].axis("off")
excluded = hours[~usable]
ax[0].vlines(excluded, 0, 130, color="lightgrey", lw=0.6, label="Excluded
records")
ax[0].plot(hours[usable], radiance[usable][:, band_window].mean(axis=1), ".",
ms=2, color="tab:blue", label="800-1000 cm$^{-1}$ window")
ax[0].plot(hours[usable], radiance[usable][:, band_co2].mean(axis=1), ".", ms=2,
color="tab:red", label="660-670 cm$^{-1}$ CO$_2$ band")
ax[0].set_ylabel(rad_label)
ax[0].legend(loc="center", bbox_to_anchor=(0.5, 0.70), ncol=3, markerscale=4,
fontsize=8)
ax[1].plot(hours[usable], obs.loc[usable, "bt_window"], ".", ms=2,
color="tab:blue", label="Window BT")
ax[1].plot(hours[usable], obs.loc[usable, "bt_co2"], ".", ms=2, color="tab:red",
label="CO$_2$ band BT")
ax[1].plot(hours, obs["air_temp"], color="k", lw=0.8, label="Surface air
temperature")
ax[1].set_ylabel("Temperature [K]")
ax[1].set_ylim(160, 285)
ax[1].legend(loc="lower center", ncol=3, markerscale=4, fontsize=8)
image = ax[2].imshow(rad_5min.values.T, aspect="auto", origin="lower", extent=[0,
24, wnum[0], wnum[-1]], vmin=0, vmax=130, cmap="viridis")
ax[2].set_ylabel(r"Wavenumber [cm$^{-1}$]")
ax[2].set_xlabel("Time [hours UTC, 1 June 2022]")
fig.colorbar(image, cax=axes[2, 1], label=rad_label)
ax[0].set_title("Temporal variation of the usable spectra (bottom: 5-minute
means)")
ax[2].set_xlim(0, 24)
fig.tight_layout()
fig.savefig("figures/fig3_radiance_time.png", dpi=200)
plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(hourly.index, hourly["present_pct"], color="lightgrey", label="Records
present")
ax.bar(hourly.index, hourly["usable_pct"], color="tab:blue", label="Usable
records")
ax.axhline(100, color="k", lw=0.8, ls="--")
ax.set_xticks(range(0, 24))
ax.set_xlabel("Hour [UTC]")
ax.set_ylabel("Percent of reference (%d records/hour)" % reference)
ax.set_title("Hourly data completeness")
ax.set_ylim(0, 110)
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig("figures/fig4_hourly_completeness.png", dpi=200)
plt.close(fig)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
bottom = np.zeros(len(hourly))
for name, color in zip(["usable", "invalid spectrum", "hatch not open",
"missing"], ["tab:blue", "tab:red", "tab:orange", "black"]):
ax[0].bar(hourly.index, hourly[name], bottom=bottom, color=color, label=name)
bottom = bottom + hourly[name].values
ax[0].set_xlabel("Hour [UTC]")
ax[0].set_ylabel("Number of records")
ax[0].set_title("Record status by hour")
ax[0].legend(fontsize=8, loc="lower right")
ax[1].hist(gap, bins=np.arange(10, 80, 2), color="tab:blue")
ax[1].set_yscale("log")
ax[1].set_xlabel("Time between consecutive records [s]")
ax[1].set_ylabel("Count")
ax[1].set_title("Time-step distribution")
fig.tight_layout()
fig.savefig("figures/fig5_missing_data_distribution.png", dpi=200)
plt.close(fig)
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
ax[0, 0].semilogy(hours, obs["anomaly"], ".", ms=2, color="tab:blue")
ax[0, 0].semilogy(hours[obs["flag_invalid"].values], obs.loc[obs["flag_invalid"],
"anomaly"], "o", ms=4, color="tab:red", label="Invalid")
ax[0, 0].axhline(suspect_limit, color="orange", ls="--", label="Suspect limit")
ax[0, 0].axhline(glitch_limit, color="red", ls="--", label="Glitch limit")
ax[0, 0].set_xlabel("Hour [UTC]")
ax[0, 0].set_ylabel("RMS departure from local median [mW/(m$^2$ sr cm$^{-1}$)]")
ax[0, 0].set_title("(a) Spectrum anomaly metric, 420-700 cm$^{-1}$")
ax[0, 0].legend(fontsize=8)
ax[0, 1].hist(d.clip(-30, 10), bins=np.arange(-30, 10.5, 0.5), color="tab:blue")
ax[0, 1].axvline(diff_median - 5 * diff_sigma, color="red", ls="--")
ax[0, 1].axvline(diff_median + 5 * diff_sigma, color="red", ls="--",
label=r"$\pm5\sigma$")
ax[0, 1].set_yscale("log")
ax[0, 1].set_xlabel("BT(660-670 cm$^{-1}$) minus air temperature [K]")
ax[0, 1].set_ylabel("Count")
ax[0, 1].set_title("(b) Physical consistency check (hatch open)")
ax[0, 1].legend(fontsize=8)
edges = [400, 600, 800, 1000, 1200, 1400, 1600, 1800]
neg_by_band = [(radiance[usable][:, (wnum >= lo) & (wnum < hi)] < 0).mean() * 100
for lo, hi in zip(edges[:-1], edges[1:])]
ax[1, 0].bar(range(len(neg_by_band)), neg_by_band, color="tab:blue")
ax[1, 0].set_xticks(range(len(neg_by_band)))
ax[1, 0].set_xticklabels(["%d-%d" % (lo, hi) for lo, hi in zip(edges[:-1],
edges[1:])], rotation=30)
ax[1, 0].set_xlabel(r"Spectral band [cm$^{-1}$]")
ax[1, 0].set_ylabel("Negative values (%)")
ax[1, 0].set_title("(c) Negative radiances in usable spectra")
counts = obs["category"].value_counts().reindex(categories)
ax[1, 1].barh(counts.index, counts.values, color=["tab:blue", "tab:red",
"tab:orange", "black"])
for i, v in enumerate(counts.values):
ax[1, 1].text(v, i, " %d (%.1f%%)" % (v, 100 * v / len(obs)), va="center",
fontsize=8)
ax[1, 1].set_xlim(0, 4800)
ax[1, 1].invert_yaxis()
ax[1, 1].set_xlabel("Number of records")
ax[1, 1].set_title("(d) Summary of record status")
fig.tight_layout()
fig.savefig("figures/fig6_qc_summary.png", dpi=200)
plt.close(fig)
