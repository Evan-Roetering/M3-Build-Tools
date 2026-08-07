"""
plot_data.py

Module in membpy for plotting data from other modules
"""

def plot_raft(lookup, raft, title, inp_fig=None, raft_color="red", nonraft_color="blue"):
    """
    Make a plotly plot of a raft
    """

    import plotly.graph_objects as go

    if inp_fig is not None:
        if not isinstance(inp_fig, (go.Figure, go.FigureWidget)):
            raise TypeError("`fig` must be a plotly.graph_objects.Figure")
        fig = go.Figure(inp_fig.to_dict())
    else:
        fig = go.Figure()

    xs = []
    ys = []
    zs = []
    labels = []
    raft_xs = []
    raft_ys = []
    raft_zs = []
    raft_labels = []
    for lipid in lookup.keys():
        if lipid in raft:
            raft_xs.append(lookup[lipid]["coord"][0])
            raft_ys.append(lookup[lipid]["coord"][1])
            raft_zs.append(lookup[lipid]["coord"][2])    
            raft_labels.append(f"{lookup[lipid]['resname']}-{lipid}")
        else:
            xs.append(lookup[lipid]["coord"][0])
            ys.append(lookup[lipid]["coord"][1])
            zs.append(lookup[lipid]["coord"][2])
            labels.append(f"{lookup[lipid]['resname']}-{lipid}")

    raft_args = {"x": raft_xs, "y": raft_ys, "z": raft_zs, "name": "Raft", "mode": "markers", "marker_color": raft_color, "hovertext": raft_labels}
    nonraft_args = {"x": xs, "y": ys, "z": zs, "name": "Other Lipids", "mode": "markers", "marker_color": nonraft_color, "hovertext": labels, "opacity": 0.2}

    fig.add_trace(go.Scatter3d(**nonraft_args))
    fig.add_trace(go.Scatter3d(**raft_args))
    
    # build layout        
    scene_settings = dict(aspectmode='data')
    layout_args = {"scene"      : scene_settings,
                   "margin"     : dict(l=0, r=0, b=0, t=40),
                   "template"   : "plotly_white",
                   "showlegend" : True}
    if title is not None:
        layout_args["title"] = title
    
    # overwrite figures previous layout
    fig.update_layout(**layout_args)

    return fig

def plot_system(system, title, color="lightblue", sys_name="System", inp_fig=None):
    """
    Make a plotly plot of a raft
    """

    import plotly.graph_objects as go

    if inp_fig is not None:
        if not isinstance(inp_fig, (go.Figure, go.FigureWidget)):
            raise TypeError("`fig` must be a plotly.graph_objects.Figure")
        fig = go.Figure(inp_fig.to_dict())
    else:
        fig = go.Figure()

    xs = []
    ys = []
    zs = []
    labels = []
    for bead in system:
        xs.append(bead["coord"][0])
        ys.append(bead["coord"][1])
        zs.append(bead["coord"][2])
        labels.append(f"{bead['resname']}-{bead['resid']}")
    args = {"x": xs, "y": ys, "z": zs, "name": sys_name, "mode": "markers", "marker_color": color, "hovertext": labels}

    fig.add_trace(go.Scatter3d(**args))
    
    # build layout        
    scene_settings = dict(aspectmode='data')
    layout_args = {"scene"      : scene_settings,
                   "margin"     : dict(l=0, r=0, b=0, t=40),
                   "template"   : "plotly_white",
                   "showlegend" : True}
    if title is not None:
        layout_args["title"] = title
    
    # overwrite figures previous layout
    fig.update_layout(**layout_args)

    return fig

def plot_positions(positions, title, color="lightblue", sys_name="System", inp_fig=None):
    """
    Make a plotly plot of a raft
    """

    import plotly.graph_objects as go

    if inp_fig is not None:
        if not isinstance(inp_fig, (go.Figure, go.FigureWidget)):
            raise TypeError("`fig` must be a plotly.graph_objects.Figure")
        fig = go.Figure(inp_fig.to_dict())
    else:
        fig = go.Figure()

    args = {"x": positions[:,0], "y": positions[:,1], "z": positions[:,2], "name": sys_name, "mode": "markers", "marker_color": color}

    fig.add_trace(go.Scatter3d(**args))
    
    # build layout        
    scene_settings = dict(aspectmode='data')
    layout_args = {"scene"      : scene_settings,
                   "margin"     : dict(l=0, r=0, b=0, t=40),
                   "template"   : "plotly_white",
                   "showlegend" : True}
    if title is not None:
        layout_args["title"] = title
    
    # overwrite figures previous layout
    fig.update_layout(**layout_args)

    return fig

def plot_lookup(lookup, title, sys_name="System", inp_fig=None, color="red", opacity=1):
    """
    Make a plotly plot of a raft
    """

    import plotly.graph_objects as go

    if inp_fig is not None:
        if not isinstance(inp_fig, (go.Figure, go.FigureWidget)):
            raise TypeError("`fig` must be a plotly.graph_objects.Figure")
        fig = go.Figure(inp_fig.to_dict())
    else:
        fig = go.Figure()

    xs = []
    ys = []
    zs = []
    labels = []
    for lipid in lookup.keys():
        xs.append(lookup[lipid]["coord"][0])
        ys.append(lookup[lipid]["coord"][1])
        zs.append(lookup[lipid]["coord"][2])
        labels.append(f"{lookup[lipid]['resname']}-{lipid}")
    args = {"x": xs, "y": ys, "z": zs, "name": sys_name, "mode": "markers", "marker_color": color, "hovertext": labels, "opacity": opacity}

    fig.add_trace(go.Scatter3d(**args))
    
    # build layout        
    scene_settings = dict(aspectmode='data')
    layout_args = {"scene"      : scene_settings,
                   "margin"     : dict(l=0, r=0, b=0, t=40),
                   "template"   : "plotly_white",
                   "showlegend" : True}
    if title is not None:
        layout_args["title"] = title
    
    # overwrite figures previous layout
    fig.update_layout(**layout_args)

    return fig

def plot_indices(lookup, indices, title, sys_name="System", inp_fig=None, color="red", opacity=1):
    import plotly.graph_objects as go

    if inp_fig is not None:
        if not isinstance(inp_fig, (go.Figure, go.FigureWidget)):
            raise TypeError("`fig` must be a plotly.graph_objects.Figure")
        fig = go.Figure(inp_fig.to_dict())
    else:
        fig = go.Figure()

    xs = []
    ys = []
    zs = []
    labels = []
    for lipid in indices:
        xs.append(lookup[lipid]["coord"][0])
        ys.append(lookup[lipid]["coord"][1])
        zs.append(lookup[lipid]["coord"][2])
        labels.append(f"{lookup[lipid]['resname']}-{lipid}")
    args = {"x": xs, "y": ys, "z": zs, "name": sys_name, "mode": "markers", "marker_color": color, "hovertext": labels, "opacity": opacity}

    fig.add_trace(go.Scatter3d(**args))
    
    # build layout        
    scene_settings = dict(aspectmode='data')
    layout_args = {"scene"      : scene_settings,
                   "margin"     : dict(l=0, r=0, b=0, t=40),
                   "template"   : "plotly_white",
                   "showlegend" : True}
    if title is not None:
        layout_args["title"] = title
    
    # overwrite figures previous layout
    fig.update_layout(**layout_args)

    return fig

def plot_time_series(datadict, time_list, title, xlabel, ylabel, vline=None, vline_label=None):
    import plotly.graph_objects as go
    fig = go.Figure()
    for key in datadict.keys():
        fig.add_trace(go.Scatter(x=time_list, y=datadict[key], name=key, mode="lines", line=dict(width=4)))
    if vline is not None and vline_label is not None:
        fig.add_vline(vline, line_width=4, line_dash="dash", line_color="grey", annotation_text=vline_label, annotation_position="top right")
    elif vline is not None:
        fig.add_vline(vline, line_width=4, line_dash="dash", line_color="grey")
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel, template="plotly_white")
    # make sure that graph starts at y=0
    fig.update_yaxes(range=[0, None])
    fig.update_layout(xaxis=dict(title=dict(text=f"<b>{xlabel}</b>", font=dict(size=35, family="Arial Black"))), yaxis=dict(title=dict(text=f"<b>{ylabel}</b>", font=dict(size=35, family="Arial"))))
    fig.update_layout(legend=dict(font=dict(size=32, family="Arial")))
    fig.update_xaxes(tickfont=dict(size=22, family="Arial Black"))
    fig.update_yaxes(tickfont=dict(size=22, family="Arial Black"))
    return fig

def plot_composition(datadict, title, xlabel, ylabel):
    import plotly.graph_objects as go
    fig = go.Figure()
    for key in datadict.keys():
        fig.add_trace(go.Bar(x=[key], y=[datadict[key]], name=key))
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel, template="plotly_white")
    fig.update_layout(xaxis=dict(title=dict(text=f"<b>{xlabel}</b>", font=dict(size=35, family="Arial Black"))), yaxis=dict(title=dict(text=f"<b>{ylabel}</b>", font=dict(size=35, family="Arial"))))
    fig.update_layout(legend=dict(font=dict(size=32, family="Arial")))
    fig.update_xaxes(tickfont=dict(size=22, family="Arial Black"))
    fig.update_yaxes(tickfont=dict(size=22, family="Arial Black"))
    return fig

def interaction_matrix(interaction_dict, title, xlabel, ylabel, z_range=[-1, 1]):
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np

    residues = list(interaction_dict.keys())
    df = pd.DataFrame([[interaction_dict[r1][r2] for r2 in residues] for r1 in residues], index=residues, columns=residues)
    fig_dict = {'data': [], 'layout': {}}
    fig_dict['layout']['xaxis'] = {'title': xlabel}
    fig_dict['layout']['yaxis'] = {'title': ylabel}
    fig_dict['layout']['title'] = title
    text = np.round(df.values, 3).astype(str)
    fig_dict['data'] = go.Heatmap(x=df.columns, y=df.index, z=df.values, 
                                  zmin=z_range[0], zmax=z_range[1], 
                                  colorscale="bluered", colorbar={'title': "Log2 Interaction"}, 
                                  text=text, texttemplate="%{text}", textfont={"size":12, "color":"black"}, 
                                  hovertemplate="Interaction: %{y} → %{x}<br>Fraction: %{z:.3f}<extra></extra>")
    
    fig = go.Figure(fig_dict)

    fig.update_layout(template="plotly_white")
    fig.update_layout(xaxis=dict(title=dict(text=f"<b>{xlabel}</b>", font=dict(size=35, family="Arial Black"))), yaxis=dict(title=dict(text=f"<b>{ylabel}</b>", font=dict(size=35, family="Arial"))))
    fig.update_layout(legend=dict(font=dict(size=32, family="Arial")))
    fig.update_xaxes(tickfont=dict(size=22, family="Arial Black"))
    fig.update_yaxes(tickfont=dict(size=22, family="Arial Black"))

    return fig

def plot_rdf(rdf_dict, bin_centers, title, xlabel, ylabel):
    import plotly.graph_objects as go
    fig = go.Figure()
    for key in rdf_dict.keys():
        fig.add_trace(go.Scatter(x=bin_centers, y=rdf_dict[key], name=key, mode="lines", line=dict(width=8)))
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel, template="plotly_white")
    fig.update_layout(xaxis=dict(title=dict(text=f"<b>{xlabel}</b>", font=dict(size=35, family="Arial Black"))), yaxis=dict(title=dict(text=f"<b>{ylabel}</b>", font=dict(size=35, family="Arial"))))
    fig.update_layout(legend=dict(font=dict(size=32, family="Arial")))
    fig.update_xaxes(tickfont=dict(size=22, family="Arial Black"))
    fig.update_yaxes(tickfont=dict(size=22, family="Arial Black"))
    return fig

def box_plot(data_dict, title, xlabel, ylabel):

    import plotly.graph_objects as go
    import numpy as np

    fig = go.Figure()

    for key, values in data_dict.items():
        vals = np.asarray(values, dtype=float)

        q1 = float(np.percentile(vals, 25))
        median = float(np.percentile(vals, 50))
        q3 = float(np.percentile(vals, 75))

        iqr = q3 - q1
        lower_fence = float(np.min(vals))
        upper_fence = float(np.max(vals))

        fig.add_trace(go.Box(
            name=key,
            x=[key],
            boxpoints=False,
            q1=[q1],
            median=[median],
            q3=[q3],
            lowerfence=[lower_fence],
            upperfence=[upper_fence],
            ))

    fig.update_layout(
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        template="plotly_white",
    )
    return fig

def plot_pdf(pdf_dict, bin_count, title, xlabel, ylabel, meanline=False):
    import plotly.graph_objects as go
    import numpy as np

    # Make bins - shared across all datasets
    min_val = min([np.min(pdf_dict[key]) for key in pdf_dict.keys()])
    max_val = max([np.max(pdf_dict[key]) for key in pdf_dict.keys()])
    bins = np.linspace(min_val, max_val, bin_count)
    
    # Compute bin centers (for x-axis) and bin width (for normalization)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_width = bins[1] - bins[0]

    fig = go.Figure()
    for key in pdf_dict.keys():
        # Compute histogram counts in each bin
        counts, _ = np.histogram(pdf_dict[key], bins=bins)
        # Normalize: divide by (total count × bin_width) so PDF integrates to 1
        pdf = counts / (np.sum(counts) * bin_width)
        fig.add_trace(go.Scatter(x=bin_centers, y=pdf, name=key, mode="lines", line=dict(width=8)))
        if meanline:
            mean_val = np.mean(pdf_dict[key])
            fig.add_vline(x=mean_val, line_width=4, line_dash="dash", line_color="grey", annotation_text=f"{mean_val:.3f}", annotation_position="top right")
            
    
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel, template="plotly_white")
    fig.update_layout(xaxis=dict(title=dict(text=f"<b>{xlabel}</b>", font=dict(size=35, family="Arial Black"))), yaxis=dict(title=dict(text=f"<b>{ylabel}</b>", font=dict(size=35, family="Arial"))))
    fig.update_layout(legend=dict(font=dict(size=32, family="Arial")))
    fig.update_xaxes(tickfont=dict(size=22, family="Arial Black"))
    fig.update_yaxes(tickfont=dict(size=22, family="Arial Black"))
    return fig

def plot_lookup_bykey(lookup, title, color_key, sys_name="System", inp_fig=None, colorscale="Thermal", opacity=1, colorrange=None, show_bar=True):
    """
    Make a plotly plot of a raft with points colored by the value of a specified key in the lookup dictionary
    """

    import plotly.graph_objects as go

    if inp_fig is not None:
        if not isinstance(inp_fig, (go.Figure, go.FigureWidget)):
            raise TypeError("`fig` must be a plotly.graph_objects.Figure")
        fig = go.Figure(inp_fig.to_dict())
    else:
        fig = go.Figure()

    xs = []
    ys = []
    zs = []
    cs = []
    labels = []
    for lipid in lookup.keys():
        xs.append(lookup[lipid]["coord"][0])
        ys.append(lookup[lipid]["coord"][1])
        zs.append(lookup[lipid]["coord"][2])
        cs.append(lookup[lipid][color_key])
        labels.append(f"{lookup[lipid]['resname']}-{lipid}: {color_key}={lookup[lipid][color_key]}")
    args = {"x": xs, "y": ys, "z": zs, 
            "marker": dict(color=cs, 
                           colorscale=colorscale, 
                           cmin=colorrange[0] if colorrange else None, 
                           cmax=colorrange[1] if colorrange else None, 
                           showscale=show_bar), 
            "name": sys_name, 
            "mode": "markers", 
            "hovertext": labels, 
            "opacity": opacity}

    fig.add_trace(go.Scatter3d(**args))
    
    # build layout        
    scene_settings = dict(aspectmode='data')
    layout_args = {"scene"      : scene_settings,
                   "margin"     : dict(l=0, r=0, b=0, t=40),
                   "template"   : "plotly_white",
                   "showlegend" : True}
    if title is not None:
        layout_args["title"] = title
    
    # overwrite figures previous layout
    fig.update_layout(**layout_args)

    return fig