import numpy as np
import random
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

# Create a random list of values.

bp_dict = {"A": 1, "T": 2, "C": 3, "G": 4}


def encode_sequence(sequence):
    encoded_sequence = []
    for b in sequence:
        encoded_sequence.append(bp_dict[b])
    return encoded_sequence


sequences = [np.random.randint(1, 5, size=(1000)) for _ in range(3)]

sequence_names = ["Sequence 1", "Sequence 2", "Sequence 3"]

colorscale = ['red', 'purple', 'green', 'orange']

fig = go.Figure(data=go.Heatmap(z=sequences,
                                y=[1, 2, 3],
                                hoverongaps=False,
                                xgap=1,
                                ygap=10,
                                colorscale=colorscale))

ay = 4

fig.add_shape(
    type="rect",
    x0=300,
    y0=0,
    x1=320,
    y1=4,
    opacity=0.9,
    fillcolor="orange",
    line_color="orange",
)

fig.add_trace(
    go.Scatter(x=[555, 565, 565, 555, 555],
               y=[0, 0, ay, ay, 0],
               name="Variable Region 1"))

fig.add_trace(
    go.Scatter(x=[888, 898, 898, 888, 888],
               y=[0, 0, ay, ay, 0],
               name="Variable Region 2"))

fig.update_layout(showlegend=False)

fig.update_layout(title_text='wildtype variants')

fig.write_html("index.html")