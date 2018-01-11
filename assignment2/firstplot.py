import plotly as py
import plotly.graph_objs as go
import pandas as pd

#preparing data set for question 1
df_nutrition = pd.read_csv('Nutrition.csv')
df_death = pd.read_csv('Leading_Causes_Death.csv')
#only keep data from 2011 onwards
df_death_final = df_death.loc[df_death['Year'] >= 2011]
#remove 'All Causes'
df_death_final = df_death.loc[df_death['Cause Name'] != 'All Causes']
#remove irrelevant columns
df_death_final = df_death_final.drop('113 Cause Name', axis=1)
df_death_final = df_death_final.drop('Age-adjusted Death Rate', axis=1)
df_death_final = df_death_final.drop('Year', axis=1)
#aggregate per cause of death per state and compute average number of deaths
df_death_final = df_death_final.groupby(['Cause Name', 'State'],as_index=False).mean()

stateList = df_death_final["State"].unique()

data = []

for state in stateList:
    trace = go.Scatter(
        x=df_death_final.loc[df_death_final['State'] == state]["Cause Name"],
        y=df_death_final.loc[df_death_final['State'] == state]["Deaths"],
        name=state,
        mode='markers'
    )
    data.append(trace)


layout = go.Layout(
   title='Bla',
   xaxis=dict(title='Cause'),
   yaxis=dict(title='Avg')
   )

fig = go.Figure(data=data, layout=layout)
py.offline.plot(fig)