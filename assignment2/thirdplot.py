import plotly as py
import plotly.graph_objs as go
import pandas as pd
import numpy as np

#preparing data set for question 3

df = pd.read_csv('Nutrition.csv')

#only keep interesting columns
df_nutrition = df[['YearStart','LocationDesc','Education','Income','Gender','Age(years)','Race/Ethnicity']]

#assume clicked state
df_nutrition = df_nutrition.loc[df_nutrition['LocationDesc'] == 'Alabama']

#here we need something like dropdown list
genderList = df_nutrition["Gender"].unique()
ageList = df_nutrition["Age(years)"].unique()

df_nutrition = df_nutrition.replace(np.nan, 'unknown')

#assume the user selected from drop down list and keep selected
#df_nutrition = df_nutrition.loc[df_nutrition['Gender'] == 'unknown']#'Male']
#df_nutrition = df_nutrition.loc[df_nutrition['Age(years)'] == 'unknown']#'18 - 24']

#prepare dataset income
df_nutrition.replace('Less than $15,000', '25k')
df_nutrition.replace('$15,000 - $24,999', '25k')
df_nutrition.replace('$25,000 - $34,999', '$25k-$50k')
df_nutrition.replace('$35,000 - $49,999', '$25k-$50k')
df_nutrition.replace('$50,000 - $74,999', '>$50k')
df_nutrition.replace('$75,000 or greater', '>$50k')

#convert categories into numerical vals
df_nutrition['Education'] = pd.Categorical(df_nutrition['Education'])
df_nutrition['EducationCode'] = df_nutrition['Education'].cat.codes
df_nutrition['Gender'] = pd.Categorical(df_nutrition['Gender'])
df_nutrition['GenderCode'] = df_nutrition['Gender'].cat.codes
df_nutrition['Income'] = pd.Categorical(df_nutrition['Income'])
df_nutrition['IncomeCode'] = df_nutrition['Income'].cat.codes
df_nutrition['Age(years)'] = pd.Categorical(df_nutrition['Age(years)'])
df_nutrition['Age(years)Code'] = df_nutrition['Age(years)'].cat.codes
df_nutrition['Race/Ethnicity'] = pd.Categorical(df_nutrition['Race/Ethnicity'])
df_nutrition['Race/EthnicityCode'] = df_nutrition['Race/Ethnicity'].cat.codes

# data = [
#     go.Parcoords(
#         line = dict(color = df_nutrition['YearStart'],
#                     colorscale=[[0, '#D7C16B'], [0.5, '#23D8C3'], [1, '#F3F10F']]),
#         dimensions = list([
#             dict(range = [0,800],
#                 #constraintrange = [4,8],
#                 label = 'Less than highschool', values = df_nutrition.loc[df_nutrition['Education'] == 'Less than high school'].count()),
#             dict(range = [0,8000],
#                 label = 'Highschool', values = df_nutrition.loc[df_nutrition['Education'] == 'High school graduate'].count()),
#             dict(range = [0,800],
#                 label = 'Technical school', values = df_nutrition.loc[df_nutrition['Education'] == 'Some college or technical school'].count()),
#             dict(range = [0,8000],
#                 label = 'College', values = df_nutrition.loc[df_nutrition['Education'] == 'College graduate'].count()),
#             dict(range = [0,800],
#                 label = '<$25k', values = df_nutrition.loc[df_nutrition['Income'] == '<$25k'].count()),
#             dict(range = [0,800],
#                 label = '$25k-$50k', values = df_nutrition.loc[df_nutrition['Education'] == '$25k-$50k'].count()),
#             dict(range = [0,800],
#                 label = '>$50k', values = df_nutrition.loc[df_nutrition['Education'] == '>$50k'].count())
#         ])
#     )
# ]

data = [
    go.Parcoords(
        line = dict(color = df_nutrition['YearStart'],
                    colorscale=[[0, '#D7C16B'], [0.5, '#23D8C3'], [1, '#F3F10F']]),
        dimensions = list([
            dict(range = [0,10],
                label = 'school', values = df_nutrition['EducationCode']),
            dict(range = [0,10],
                label = 'income', values = df_nutrition['IncomeCode']),
            dict(range = [0,10],
                label = 'gender', values = df_nutrition['GenderCode']),
            dict(range = [0,10],
                label = 'age', values = df_nutrition['Age(years)Code']),
            dict(range = [0,10],
                label = 'race', values = df_nutrition['Race/EthnicityCode'])
        ])
    )
]

layout = go.Layout(
    title='Bla',
    plot_bgcolor = '#E5E5E5',
    paper_bgcolor = '#E5E5E5'
)

fig = go.Figure(data=data, layout=layout)
py.offline.plot(fig)