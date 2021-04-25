import json
import numpy as np
import random
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

from sequence import sc

# bp_dict = {"A": 1, "T": 2, "C": 3, "G": 4}

# title_text = "wildtype variants"

# sequences = [
#     np.random.choice([1., 2., 3., 4., 5.], size=(2200))
#     for _ in range(num_sequences)
# ]

# sequences[0][50:55] = np.nan

# def encode_sequence(sequence):
#     encoded_sequence = []
#     for b in sequence:
#         encoded_sequence.append(bp_dict[b])
#     return encoded_sequence


def make_plotly(sequences, reference_sequence, annotation_file, reference_name,
                title, output_file):
    '''
        List of numpy arrays in 1,2,3,4,NaN form
        Last sequence is the reference
    
    '''

    # compute mutations

    ref_stack = np.vstack(
        [reference_sequence for _ in range(sequences.shape[0])])
    print(ref_stack.shape)
    sequences[np.where(sequences == ref_stack)] = np.nan

    num_sequences = len(sequences) + 1

    sequences = np.append(sequences, reference_sequence.reshape(1, -1), axis=0)

    sequence_names = [f"Sequence {str(i)}"
                      for i in range(num_sequences - 1)] + [reference_name]

    colorscale = ['grey', 'red', 'purple', 'green', 'orange', 'grey', 'grey']

    sequence_ticks = list(np.linspace(1, num_sequences, num_sequences))

    fig = go.Figure(data=go.Heatmap(z=sequences,
                                    y=sequence_ticks,
                                    hoverongaps=False,
                                    xgap=1,
                                    ygap=1,
                                    colorscale=colorscale,
                                    showscale=False))

    fig.update_layout(yaxis=dict(
        tickmode='array', tickvals=sequence_ticks, ticktext=sequence_names))

    with open(annotation_file) as fp:
        # a list of annotations
        annotations = json.load(fp)["annotations"]

    for annotation in annotations:
        x0 = annotation["DNA_start"]
        x1 = annotation["DNA_end"]
        y0 = 0.4
        if "VR " in annotation["name"]:
            y1 = num_sequences + 0.6
            color = "darkblue"
        else:
            y1 = num_sequences + 1.
            color = "darkviolet"
        fig.add_shape(type="rect",
                      x0=x0,
                      y0=y0,
                      x1=x1,
                      y1=y1,
                      fillcolor=color,
                      line={"width": 0},
                      layer="below",
                      opacity=.2)
        fig.add_trace(
            go.Scatter(x=[(x0 + x1) / 2],
                       y=[y1],
                       text=annotation["name"],
                       name=annotation["notes"],
                       hoverinfo="name",
                       hoverlabel={"namelength": -1},
                       mode="text",
                       textposition="top center"))

    fig.update_xaxes(title="Nucleotide", range=[0, 1000])
    fig.update_layout(showlegend=False)
    fig.update_layout(title_text=title)
    fig.write_html(output_file)


annotation_file = "annotations/AAV_capsid.json"
output_file = "annotations.html"

make_plotly(sc.num_stack,
            sc.reference.num,
            annotation_file=annotation_file,
            reference_name="AAV6",
            title="prelim capsids",
            output_file=output_file)
