import plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import Scatter, Layout
import pandas as pd

#plotly.offline.plot({"data": [Scatter(x=[1,2,3,4],y=[4,3,2,1])],"layout": Layout(title="hello")})

df = pd.read_csv('Nutrition.csv')
#print(df)

data = [go.Histogram(x=df["LocationAbbr"])]

layout = go.Layout(
    title='Nutrition',
    xaxis=dict(title='Location'),
    yaxis=dict(title='Count')
    )

fig = go.Figure(data=data, layout=layout)
py.offline.plot(fig)