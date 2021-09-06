from collections import defaultdict
from ..utils.definitions import *
from ..utils.utils import *
import plotly.graph_objects as go


def make_snp_plots(msdata, markers, out_prefix, data_name):
    """
    Writes SNP quality plots for the markers
    """
    markers_list = list(markers.keys())
    for marker in markers_list:
        if marker not in msdata:
            continue
        marker_data = msdata[marker].data
        samples = list(marker_data.keys())
        coords = defaultdict(dict)
        for sample in samples:
            allele = marker_data[sample].__str__()
            xvalue = marker_data[sample].get_xvalue()
            yvalue = marker_data[sample].get_yvalue()
            if allele not in coords:
                coords[allele]['x'] = []
                coords[allele]['y'] = []
            coords[allele]['x'].append(xvalue)
            coords[allele]['y'].append(yvalue)
        alleles = [markers[marker].get_allele_x(
        ), markers[marker].get_allele_y(), markers[marker].get_allele_xy()]
        alleles += MISSING_ALLELES
        out_image = out_prefix + '_' + data_name + '_' + marker + ".png"
        fig = go.Figure()
        fig.update_layout(
            title={
                'text': data_name + " - " + marker,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top',
                'font': dict(
                    family="Courier New, monospace",
                    size=24,
                    color="RebeccaPurple")
            },
            height=600,
            width=800,
        )
        i = 0
        xLimit = 0
        yLimit = 0
        for allele in alleles:
            if allele not in coords:
                continue
            xList = list(map(num, coords[allele]['x']))
            yList = list(map(num, coords[allele]['y']))
            xMax = max(xList)
            yMax = max(yList)
            if xMax > xLimit:
                xLimit = xMax + 0.5
            if yMax > yLimit:
                yLimit = yMax + 0.5
            # xAverage = round(average(xList), 4)
            # yAverage = round(average(yList), 4)
            fig.add_trace(go.Scatter(
                x=xList,
                y=yList,
                mode="markers",
                name=allele,
                opacity=0.6,
                marker=dict(color=COLORS[i])
            ))
            # fig.add_trace(go.Scatter(
            #     x=[xAverage],
            #     y=[yAverage],
            #     mode="markers",
            #     marker_symbol='diamond',
            #     name='Centroid - ' + allele,
            #     opacity=0.9,
            #     marker=dict(color=COLORS[i], size=10)
            # ))
            i += 1
        fig.update_xaxes(title_text='X')
        fig.update_yaxes(title_text='Y')
        fig.update_layout(xaxis_range=[0, xLimit], yaxis_range=[0, yLimit])
        fig.write_image(out_image)
