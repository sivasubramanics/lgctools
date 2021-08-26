# from wordcloud import WordCloud, STOPWORDS
import matplotlib.pyplot as plt
import pandas as pd
import sys
from classes.Getoptions import *


options = Getoptions(sys.argv)
print(sys.argv[3])


if False:
    filename = sys.argv[1]
    outfile = filename + '.png'
    df = pd.read_csv(filename, sep="\t", index_col='marker_count')
    text = " ".join(cat.split()[0] for cat in df.marker)
    stopwords = set(STOPWORDS)
    wordcloud = WordCloud(width = 800, height = 800,
                    background_color ='white',
                    stopwords = stopwords,
                    min_font_size = 10).generate(text)
    wordcloud.to_file(outfile)
