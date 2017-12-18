from xmatch import xmatch
import pandas

df = pandas.read_csv('data/table.csv', sep=';')
df.reset_index(inplace=True)

df['snr'] = df['nufnu_5keV(erg.s-1.cm-2)']/df['nufnu_error_5keV(erg.s-1.cm-2)']

cols = {'ra':'RA','dec':'DEC','id':'index'}
#cnn = xmatch(df,df,columns_A=cols, columns_B=cols, method='nn')
cgc = xmatch(df,df,columns_A=cols, columns_B=cols, radius=6/3600, method='gc', snr_column='snr')
df = df.iloc[cgc.index] 

