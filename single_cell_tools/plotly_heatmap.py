#!/usr/bin/python

import plotly
import plotly.plotly as py
from plotly.graph_objs import *
import plotly.figure_factory as FF

import numpy as np
from scipy.spatial.distance import pdist, squareform
import ipdb

def main():
	# get data
	data = np.genfromtxt("http://files.figshare.com/2133304/ExpRawData_E_TABM_84_A_AFFY_44.tab",
						 names=True,usecols=tuple(range(1,30)),dtype=float, delimiter="\t")
	data_array = data.view((np.float, len(data.dtype.names)))
	data_array = data_array.transpose()
	labels = data.dtype.names

	# Initialize figure by creating upper dendrogram
	figure = FF.create_dendrogram(data_array, orientation='bottom', labels=labels)
	for i in range(len(figure['data'])):
		figure['data'][i]['yaxis'] = 'y2'

	# Create Side Dendrogram
	dendro_side = FF.create_dendrogram(data_array, orientation='right')
	for i in range(len(dendro_side['data'])):
		dendro_side['data'][i]['xaxis'] = 'x2'

	#convert scipy dend to plotly
	 {'hoverinfo': 'text',
  'marker': {'color': 'rgb(0,116,217)'},
  'mode': 'lines',
  'text': None,
  'type': 'scatter',
  'x': array([-202.84413838, -332.62229802, -332.62229802, -249.43783233]),
  'xaxis': 'x2',
  'y': array([ 74.4921875,  74.4921875, 250.       , 250.       ]),
  'yaxis': 'y'}

	# Add Side Dendrogram Data to Figure
	figure['data'].extend(dendro_side['data'])

	# Create Heatmap
	dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
	dendro_leaves = list(map(int, dendro_leaves))
	data_dist = pdist(data_array)
	heat_data = squareform(data_dist)
	heat_data = heat_data[dendro_leaves,:]
	heat_data = heat_data[:,dendro_leaves]

	heatmap = Data([
		Heatmap(
			x = dendro_leaves,
			y = dendro_leaves,
			z = heat_data,
			colorscale = 'YIGnBu'
		)
	])

	heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
	heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

	# Add Heatmap Data to Figure
	figure['data'].extend(Data(heatmap))

	# Edit Layout
	figure['layout'].update({'width':800, 'height':800,
							 'showlegend':False, 'hovermode': 'closest',
							 })
	# Edit xaxis
	figure['layout']['xaxis'].update({'domain': [.15, 1],
									  'mirror': False,
									  'showgrid': False,
									  'showline': False,
									  'zeroline': False,
									  'ticks':""})
	# Edit xaxis2
	figure['layout'].update({'xaxis2': {'domain': [0, .15],
									   'mirror': False,
									   'showgrid': False,
									   'showline': False,
									   'zeroline': False,
									   'showticklabels': False,
									   'ticks':""}})

	# Edit yaxis
	figure['layout']['yaxis'].update({'domain': [0, .85],
									  'mirror': False,
									  'showgrid': False,
									  'showline': False,
									  'zeroline': False,
									  'showticklabels': False,
									  'ticks': ""})
	# Edit yaxis2
	figure['layout'].update({'yaxis2':{'domain':[.825, .975],
									   'mirror': False,
									   'showgrid': False,
									   'showline': False,
									   'zeroline': False,
									   'showticklabels': False,
									   'ticks':""}})

	# Plot!
	url = plotly.offline.plot(figure, filename='dendrogram_with_heatmap', validate=False, auto_open=False)
	print(url)
	#~ py.iplot(figure, )
	
def plot_tree( P, sel_color=None ):
    ipdb.set_trace()
    icoord = sc.array( P['icoord'] )
    dcoord = sc.array( P['dcoord'] )
    color_list = sc.array( P['color_list'] )
    xmin, xmax = icoord.min(), icoord.max()
    ymin, ymax = dcoord.min(), dcoord.max()
    itemindex = np.where(color_list==sel_color)
    
    if sel_color:
        icoord = icoord[itemindex]
        dcoord = dcoord[itemindex]
        color_list = color_list[itemindex]
    


if __name__ == "__main__":
	main()
