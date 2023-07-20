#!/usr/bin/env python3

# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
from cooltools import insulation

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.4'):
	raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
						 "please check your cooltools version and update to the latest")

import cooltools
import itertools
from matplotlib.ticker import EngFormatter

from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
from skimage.filters import threshold_li, threshold_otsu
import bbi

def getdata():
	data_dir = './data/'
	cool_file = cooltools.download_data("HFF_MicroC", cache=True, data_dir=data_dir)

	resolution = 10000
	clr = cooler.Cooler(f'{data_dir}test.mcool::resolutions/{resolution}')
	windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
	insulation_table = insulation(clr, windows, verbose=True)

	norm = LogNorm(vmax=0.1, vmin=0.001)
	return data_dir, resolution, clr, windows, insulation_table, norm

def print_first_summary(windows, insulation_table):
	first_window_summary =insulation_table.columns[[ str(windows[-1]) in i for i in insulation_table.columns]]
	print(insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005])

# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
	start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
	n = matrix_c.shape[0]
	t = np.array([[1, 0.5], [-1, 0.5]])
	matrix_a = np.dot(np.array([
		(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)
	]), t)
	x = matrix_a[:, 1].reshape(n + 1, n + 1)
	y = matrix_a[:, 0].reshape(n + 1, n + 1)
	im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
	im.set_rasterized(True)
	return im

def format_ticks(ax, x=True, y=True, rotate=True):
	bp_formatter = EngFormatter('b')
	if y:
		ax.yaxis.set_major_formatter(bp_formatter)
	if x:
		ax.xaxis.set_major_formatter(bp_formatter)
		ax.xaxis.tick_bottom()
	if rotate:
		ax.tick_params(axis='x',rotation=45)

def data_from_clr(clr, region):
	return clr.matrix(balance=True).fetch(region)

def makeregion(chr, start, windows):
	end = start + 90 * windows[0]
	return (chr, start, end)

def plot_contacts_simple(windows, resolution, data, insulation_table, region, norm):
	plt.rcParams['font.size'] = 12

	f, ax = plt.subplots(figsize=(18, 6))
	im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
	ax.set_aspect(0.5)
	ax.set_ylim(0, 10*windows[0])
	format_ticks(ax, rotate=False)
	ax.xaxis.set_visible(False)

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
	plt.colorbar(im, cax=cax)

	insul_region = bioframe.select(insulation_table, region)

	ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
	ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))
	ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
				insul_region['log2_insulation_score_'+str(windows[0])],
				label=f'Window {windows[0]} bp')

	ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

	format_ticks(ins_ax, y=False, rotate=False)
	ax.set_xlim(region[1], region[2])
	plt.savefig("fig1.png")

	for res in windows[1:]:
		ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'], label=f'Window {res} bp')
	ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
	# f
	plt.savefig("fig2.png")

def plot_with_boundaries(insulation_table, region, data, resolution, norm, windows):
	f, ax = plt.subplots(figsize=(20, 10))
	im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
	ax.set_aspect(0.5)
	ax.set_ylim(0, 10*windows[0])
	format_ticks(ax, rotate=False)
	ax.xaxis.set_visible(False)

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
	plt.colorbar(im, cax=cax)

	insul_region = bioframe.select(insulation_table, region)

	ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)

	ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
				insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

	boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]
	weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
	strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]
	ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
				weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
	ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
				strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

	ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);

	format_ticks(ins_ax, y=False, rotate=False)
	ax.set_xlim(region[1], region[2])
	plt.savefig("chr_with_boundaries.png")
	return boundaries, weak_boundaries, strong_boundaries

def make_histkwargs():
	histkwargs = dict(
		bins=10**np.linspace(-4,1,200),
		histtype='step',
		lw=2,
	)
	return histkwargs

def call_boundaries(insulation_table, windows, histkwargs):
	f, axs = plt.subplots(len(windows), 1, sharex=True, figsize=(6,6), constrained_layout=True)
	thresholds_li = {}
	thresholds_otsu = {}
	for i, (w, ax) in enumerate(zip(windows, axs)):
		ax.hist(
			insulation_table[f'boundary_strength_{w}'],
			**histkwargs
		)
		thresholds_li[w] = threshold_li(insulation_table[f'boundary_strength_{w}'].dropna().values)
		thresholds_otsu[w] = threshold_otsu(insulation_table[f'boundary_strength_{w}'].dropna().values)
		n_boundaries_li = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_li[w]).sum()
		n_boundaries_otsu = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_otsu[w]).sum()
		ax.axvline(thresholds_li[w], c='green')
		ax.axvline(thresholds_otsu[w], c='magenta')
		ax.text(0.01, 0.9,
				 f'Window {w//1000}kb',
				 ha='left',
				 va='top',
				 transform=ax.transAxes)
		ax.text(0.01, 0.7,
				f'{n_boundaries_otsu} boundaries (Otsu)',
				c='magenta',
				ha='left',
				va='top',
				transform=ax.transAxes)
		ax.text(0.01, 0.5,
				f'{n_boundaries_li} boundaries (Li)',
				c='green',
				ha='left',
				va='top',
				transform=ax.transAxes)

		ax.set(
			xscale='log',
			ylabel='# boundaries'
		)

	axs[-1].set(xlabel='Boundary strength')
	plt.savefig("boundaries.png")
	return thresholds_li, thresholds_otsu, n_boundaries_li, n_boundaries_otsu

def tabulate_boundaries(insulation_table, windows):
	is_boundary = np.any([
			~insulation_table[f'boundary_strength_{w}'].isnull()
			for w in windows],
		axis=0)
	boundaries = insulation_table[is_boundary]
	# boundaries.head()
	return boundaries

def write_thresh(path, threshes):
	with open(path, "w") as fp:
		fp.write(str(threshes))

def write_nbound(path, nb):
	with open(path, "w") as fp:
		fp.write(str(nb))

def write_bound(path, b):
	b.to_csv(path, sep = '\t', index = False)
	# with open(path, "w") as fp:
	# 	fp.write(str(b))

def main():
	data_dir, resolution, clr, windows, insulation_table, norm = getdata()
	print_first_summary(windows, insulation_table)

	region = makeregion("chr2", 10_500_000, windows)
	data = data_from_clr(clr, region)

	boundaries, weak_boundaries, strong_boundaries = plot_with_boundaries(insulation_table, region, data, resolution, norm, windows)
	write_bound("boundaries.txt", boundaries)
	write_bound("weak_boundaries.txt", weak_boundaries)
	write_bound("strong_boundaries.txt", strong_boundaries)

	thresholds_li, thresholds_otsu, n_boundaries_li, n_boundaries_otsu = call_boundaries(insulation_table, windows, make_histkwargs())
	tabulate_boundaries(insulation_table, windows)
	write_thresh("thresholds_li.txt", thresholds_li)
	write_thresh("thresholds_otsu.txt", thresholds_otsu)
	write_nbound("n_boundaries_li.txt", n_boundaries_li)
	write_nbound("n_boundaries_otsu.txt", n_boundaries_otsu)

if __name__ == "__main__":
	main()
