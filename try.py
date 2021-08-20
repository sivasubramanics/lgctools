import pandas as pd

filename = "output/BestMarkerSummaryInd.txt"
df = pd.read_csv(filename, sep="\t", index_col='marker_count')
count = 53
print(list(df.loc[count]['marker']))

