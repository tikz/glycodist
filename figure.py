from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


d = pd.read_csv("result.csv")
df = d.groupby(["UniProtID", "Pos"]).mean()
f, ax = plt.subplots(figsize=(8, 8))
ax = sns.distplot(df["ClosestGlycoDistance"], kde_kws={"bw": 0.2})
ax = sns.distplot(df["2ndClosestGlycoDistance"], kde_kws={"bw": 0.2})
ax = sns.distplot(df["FurthestResDistance"], kde_kws={"bw": 0.2})
ax.set(xlabel="Distance (Å)", ylabel="Probability density")
ax.set_xlim(0, 180)
f.savefig("figure.pdf")

print(d.groupby(["UniProtID"]).mean().shape)
print(d.groupby(["PDB ID"]).mean().shape)
percentile = stats.percentileofscore(df["ClosestGlycoDistance"], 70)
print(percentile)
