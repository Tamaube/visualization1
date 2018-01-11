import plotly as py
import plotly.graph_objs as go
import pandas as pd

#preparing data set for question 2

df = pd.read_csv('Nutrition.csv')

#only keep interesting columns
df_nutrition = df[['YearStart','LocationDesc','Question','Data_Value','Gender']]

#assume clicked state
df_nutrition = df_nutrition.loc[df_nutrition['LocationDesc'] == 'Alabama']

#here we need something like dropdown list
questionList = df_nutrition["Question"].unique()

#selected from list
df_nutrition = df_nutrition.loc[df_nutrition['Question'] == 'Percent of adults aged 18 years and older who have obesity']

data = []

genderList = df_nutrition["Gender"].unique()

for gender in genderList:
    trace = go.Bar(
        x=df_nutrition.loc[df_nutrition['Gender'] == gender]["YearStart"],
        y=df_nutrition.loc[df_nutrition['Gender'] == gender]["Data_Value"],
        name=gender
    )
    data.append(trace)

layout = go.Layout(
    barmode='stack',
    title='Bla',
    xaxis=dict(title='year'),
    yaxis=dict(title='percentage')
)

fig = go.Figure(data=data, layout=layout)
py.offline.plot(fig)