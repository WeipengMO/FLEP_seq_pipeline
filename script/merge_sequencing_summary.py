import pandas as pd
import glob


outfn = 'sequencing_summary.txt'
ss_files = [fn for fn in glob.glob('guppy_out*/*/sequencing_summary.txt')]
df_list = []
for fn in ss_files:
    df_list.append(pd.read_csv(fn, sep='\t'))

df = pd.concat(df_list, ignore_index=True)
df.to_csv(outfn, sep='\t', index=False)
