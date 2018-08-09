

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import sem

clinicaltrial_df = pd.read_csv("raw_data/clinicaltrial_data.csv")
mouse_drug_df = pd.read_csv("raw_data/mouse_drug_data.csv")
df = pd.merge(clinicaltrial_df, mouse_drug_df,how = 'left', on="Mouse ID")
df1 = df[(df['Drug']  == "Capomulin")|(df['Drug']  == 'Infubinol')|(df['Drug']  == 'Ketapril') |(df['Drug']  == 'Placebo')]
```


```python
df_t = df1.pivot_table(index='Timepoint',columns='Drug',aggfunc={'Tumor Volume (mm3)':['mean','sem']})
df_t.columns = df_t.columns.droplevel(0)
df_t
#df_t2 = df1.groupby(['Drug', 'Timepoint']).agg({'Tumor Volume (mm3)': ['mean','sem']}).unstack(level='Drug')
#df_t2.columns = g.columns.droplevel(0)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead tr th {
        text-align: left;
    }

    .dataframe thead tr:last-of-type th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th colspan="4" halign="left">mean</th>
      <th colspan="4" halign="left">sem</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>45.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>44.266086</td>
      <td>47.062001</td>
      <td>47.389175</td>
      <td>47.125589</td>
      <td>0.448593</td>
      <td>0.235102</td>
      <td>0.264819</td>
      <td>0.218091</td>
    </tr>
    <tr>
      <th>10</th>
      <td>43.084291</td>
      <td>49.403909</td>
      <td>49.582269</td>
      <td>49.423329</td>
      <td>0.702684</td>
      <td>0.282346</td>
      <td>0.357421</td>
      <td>0.402064</td>
    </tr>
    <tr>
      <th>15</th>
      <td>42.064317</td>
      <td>51.296397</td>
      <td>52.399974</td>
      <td>51.359742</td>
      <td>0.838617</td>
      <td>0.357705</td>
      <td>0.580268</td>
      <td>0.614461</td>
    </tr>
    <tr>
      <th>20</th>
      <td>40.716325</td>
      <td>53.197691</td>
      <td>54.920935</td>
      <td>54.364417</td>
      <td>0.909731</td>
      <td>0.476210</td>
      <td>0.726484</td>
      <td>0.839609</td>
    </tr>
    <tr>
      <th>25</th>
      <td>39.939528</td>
      <td>55.715252</td>
      <td>57.678982</td>
      <td>57.482574</td>
      <td>0.881642</td>
      <td>0.550315</td>
      <td>0.755413</td>
      <td>1.034872</td>
    </tr>
    <tr>
      <th>30</th>
      <td>38.769339</td>
      <td>58.299397</td>
      <td>60.994507</td>
      <td>59.809063</td>
      <td>0.934460</td>
      <td>0.631061</td>
      <td>0.934121</td>
      <td>1.218231</td>
    </tr>
    <tr>
      <th>35</th>
      <td>37.816839</td>
      <td>60.742461</td>
      <td>63.371686</td>
      <td>62.420615</td>
      <td>1.052241</td>
      <td>0.984155</td>
      <td>1.127867</td>
      <td>1.287481</td>
    </tr>
    <tr>
      <th>40</th>
      <td>36.958001</td>
      <td>63.162824</td>
      <td>66.068580</td>
      <td>65.052675</td>
      <td>1.223608</td>
      <td>1.055220</td>
      <td>1.158449</td>
      <td>1.370634</td>
    </tr>
    <tr>
      <th>45</th>
      <td>36.236114</td>
      <td>65.755562</td>
      <td>70.662958</td>
      <td>68.084082</td>
      <td>1.223977</td>
      <td>1.144427</td>
      <td>1.453186</td>
      <td>1.351726</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax = plt.subplots(figsize=(9,6))
markers = ['o','s','p','*']
lines = [':','--','-.',':']
for i in range(0,int(len(df_t.columns)/2)):
    x = df_t.index
    y = df_t.iloc[:, i]
    error = df_t.iloc[:, i+4] 
    ax.errorbar(x, y,yerr=error,ls=lines[i],lw=2,marker=markers[i], capsize=4, markersize=8,label=(list(df_t)[i])[1])

plt.ylim(30, None)
plt.xlim(0, 45)
plt.title('Tumor Response to Treatment',fontsize=18)
ax.set_ylabel(r"Tumor Volume (mm3)", fontsize=18, color="black")
ax.set_xlabel(r"Timepoint",fontsize=18, color="black")
plt.grid(linestyle = '--', linewidth = 0.3)
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')
ax.legend(loc='upper left', frameon=True)
plt.show()
```


![png](output_2_0.png)



```python
df_m = df1.pivot_table(index='Timepoint',columns='Drug',aggfunc={'Metastatic Sites':['mean','sem']})
df_m.columns = df_m.columns.droplevel(0)
df_m
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead tr th {
        text-align: left;
    }

    .dataframe thead tr:last-of-type th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th colspan="4" halign="left">mean</th>
      <th colspan="4" halign="left">sem</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.160000</td>
      <td>0.280000</td>
      <td>0.304348</td>
      <td>0.375000</td>
      <td>0.074833</td>
      <td>0.091652</td>
      <td>0.098100</td>
      <td>0.100947</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.320000</td>
      <td>0.666667</td>
      <td>0.590909</td>
      <td>0.833333</td>
      <td>0.125433</td>
      <td>0.159364</td>
      <td>0.142018</td>
      <td>0.115261</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.375000</td>
      <td>0.904762</td>
      <td>0.842105</td>
      <td>1.250000</td>
      <td>0.132048</td>
      <td>0.194015</td>
      <td>0.191381</td>
      <td>0.190221</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.652174</td>
      <td>1.050000</td>
      <td>1.210526</td>
      <td>1.526316</td>
      <td>0.161621</td>
      <td>0.234801</td>
      <td>0.236680</td>
      <td>0.234064</td>
    </tr>
    <tr>
      <th>25</th>
      <td>0.818182</td>
      <td>1.277778</td>
      <td>1.631579</td>
      <td>1.941176</td>
      <td>0.181818</td>
      <td>0.265753</td>
      <td>0.288275</td>
      <td>0.263888</td>
    </tr>
    <tr>
      <th>30</th>
      <td>1.090909</td>
      <td>1.588235</td>
      <td>2.055556</td>
      <td>2.266667</td>
      <td>0.172944</td>
      <td>0.227823</td>
      <td>0.347467</td>
      <td>0.300264</td>
    </tr>
    <tr>
      <th>35</th>
      <td>1.181818</td>
      <td>1.666667</td>
      <td>2.294118</td>
      <td>2.642857</td>
      <td>0.169496</td>
      <td>0.224733</td>
      <td>0.361418</td>
      <td>0.341412</td>
    </tr>
    <tr>
      <th>40</th>
      <td>1.380952</td>
      <td>2.100000</td>
      <td>2.733333</td>
      <td>3.166667</td>
      <td>0.175610</td>
      <td>0.314466</td>
      <td>0.315725</td>
      <td>0.297294</td>
    </tr>
    <tr>
      <th>45</th>
      <td>1.476190</td>
      <td>2.111111</td>
      <td>3.363636</td>
      <td>3.272727</td>
      <td>0.202591</td>
      <td>0.309320</td>
      <td>0.278722</td>
      <td>0.304240</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax = plt.subplots(figsize=(9,6))
markers = ['o','s','p','*']
lines = [':','--','-.',':']
for i in range(0,int(len(df_m.columns)/2)):
    x = df_m.index
    y = df_m.iloc[:, i]
    error = df_m.iloc[:, i+4] 
    ax.errorbar(x, y,yerr=error,ls=lines[i],lw=2,marker=markers[i], capsize=4, markersize=8,label=(list(df_m)[i])[1])

plt.ylim(0, None)
plt.xlim(0, 45)
plt.title('Metastatic Spread During Treatment',fontsize=18)
ax.set_ylabel(r"Met. Sites", fontsize=18, color="black")
ax.set_xlabel(r"Timepoint",fontsize=18, color="black")
plt.grid(linestyle = '--', linewidth = 0.3)
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')
ax.legend(loc='upper left', frameon=True)
plt.show()
```


![png](output_4_0.png)



```python
df_n =df1.pivot_table(index='Timepoint',columns='Drug',values='Mouse ID',aggfunc=[lambda x: x.count()*100/25])
df_n.columns = df_n.columns.droplevel(0)
df_n
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>100.0</td>
      <td>100.0</td>
      <td>100.0</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>100.0</td>
      <td>100.0</td>
      <td>92.0</td>
      <td>96.0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>100.0</td>
      <td>84.0</td>
      <td>88.0</td>
      <td>96.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>96.0</td>
      <td>84.0</td>
      <td>76.0</td>
      <td>80.0</td>
    </tr>
    <tr>
      <th>20</th>
      <td>92.0</td>
      <td>80.0</td>
      <td>76.0</td>
      <td>76.0</td>
    </tr>
    <tr>
      <th>25</th>
      <td>88.0</td>
      <td>72.0</td>
      <td>76.0</td>
      <td>68.0</td>
    </tr>
    <tr>
      <th>30</th>
      <td>88.0</td>
      <td>68.0</td>
      <td>72.0</td>
      <td>60.0</td>
    </tr>
    <tr>
      <th>35</th>
      <td>88.0</td>
      <td>48.0</td>
      <td>68.0</td>
      <td>56.0</td>
    </tr>
    <tr>
      <th>40</th>
      <td>84.0</td>
      <td>40.0</td>
      <td>60.0</td>
      <td>48.0</td>
    </tr>
    <tr>
      <th>45</th>
      <td>84.0</td>
      <td>36.0</td>
      <td>44.0</td>
      <td>44.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax = plt.subplots(figsize=(9,6))
markers = ['o','s','p','*']
lines = [':','--','-.',':']
for i in range(0,len(df_n.columns)):
    x = df_n.index
    y = df_n.iloc[:, i]
    ax.plot(x, y,ls=lines[i],lw=2,marker=markers[i])

plt.ylim(30, None)
plt.xlim(0, 45)
plt.title('Survival During Treatment',fontsize=18)
ax.set_ylabel(r"Survival Rate (%)", fontsize=18, color="black")
ax.set_xlabel(r"Timepoint",fontsize=18, color="black")
plt.grid(linestyle = '--', linewidth = 0.3)
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')
ax.legend(loc='lower left', frameon=True)
plt.show()
```


![png](output_6_0.png)



```python
df2 = df1[(df1['Timepoint']  == 0)|(df1['Timepoint']  == 45)]
df3 = pd.DataFrame(df2.groupby(['Drug','Timepoint'])['Tumor Volume (mm3)'].apply(lambda x:x.sum()/x.count())).unstack()
df3.columns = df3.columns.droplevel(0)
df3['total_tumor_change'] =(df3[df3.columns[1]]-df3[df3.columns[0]])*100/df3[df3.columns[0]]
df3
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Timepoint</th>
      <th>0</th>
      <th>45</th>
      <th>total_tumor_change</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Capomulin</th>
      <td>45.0</td>
      <td>36.236114</td>
      <td>-19.475303</td>
    </tr>
    <tr>
      <th>Infubinol</th>
      <td>45.0</td>
      <td>65.755562</td>
      <td>46.123472</td>
    </tr>
    <tr>
      <th>Ketapril</th>
      <td>45.0</td>
      <td>70.662958</td>
      <td>57.028795</td>
    </tr>
    <tr>
      <th>Placebo</th>
      <td>45.0</td>
      <td>68.084082</td>
      <td>51.297960</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax = plt.subplots(figsize=(8,5))
width = [0.8,0.8,0.8,0.8]
rects1 = ax.bar(df3.index,df3['total_tumor_change'],color=['green', 'red','red','red'],edgecolor='black',width=width)
plt.title('Total Tumor Change Over 45 Day Treatment',fontsize=18)
ax.set_ylabel("% Tumor Volume Change", fontsize=12)
plt.grid(linestyle = '--', linewidth = 0.3)
plt.grid(linestyle = '--', linewidth = 0.3)
plt.axhline(linewidth=1, color='black')

def autolabel(rects, ax):
    # Get y-axis height to calculate label position from.
    (y_bottom, y_top) = ax.get_ylim()
    y_height = y_top - y_bottom

    for rect in rects:
        height = rect.get_height()

        # Fraction of axis height taken up by this rectangle
        p_height = (height / y_height)

        # If we can fit the label above the column, do that;
        # otherwise, put it inside the column.
        if p_height > 0: # arbitrary; 95% looked good to me.
            label_position = 0.05*height
        else:
            label_position = height + (y_height * 0.15)

        ax.text(rect.get_x() + rect.get_width()/2., label_position,
                '%d' % int(height)+'%',
                ha='center', va='bottom')

autolabel(rects1, ax)

plt.show()
```


![png](output_8_0.png)

