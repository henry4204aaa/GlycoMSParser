import pandas as pd
df = pd.read_csv("MLgeneratetest1_score007_salvage_PLabeled.csv")
print(df["Structure"].value_counts().head(10))  # top-heavy?
print("Non-glycan count:", (df["Structure"]=="Non-glycan").sum())
