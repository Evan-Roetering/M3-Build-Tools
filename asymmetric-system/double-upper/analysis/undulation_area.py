import numpy as np
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import os
import undulation_area_helpers as helpers
import membpy

frame_path = r"./gro_snapshots/frame_{frame_num}.gro"
outpath = "./area"
first_frame = 0
last_frame = 400
frame_skip = 1
frame_spacing = 2.5
surface_beadnames = ['OOH', 'GL1', 'GL2', 'AM1', 'AM2', 'OH1', 'OH2', 'ROH']
continuity_dist = 1
coarse_res = 100
fine_res = 250

def make_outdir(outpath, dirname=None):
    """
    Create an output directory if it does not already exist.

    Parameters:
        outpath (str): The path where the output directory should be created.
        dirname (str): The name of the output directory to create.

    Returns:
        str: The full path to the created or existing output directory.
    """
    if dirname is None:
        outdir = outpath
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.join(outpath, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    return outdir

frame_list = np.arange(first_frame, last_frame + 1, frame_skip)
time_list = frame_list * frame_spacing
idx_list = np.arange(len(frame_list))

upper_areas = np.zeros(len(frame_list))
lower_areas = np.zeros(len(frame_list))
avg_areas = np.zeros(len(frame_list))
upper_apl = np.zeros(len(frame_list))
lower_apl = np.zeros(len(frame_list))
avg_apl = np.zeros(len(frame_list))

worker = partial(helpers.process_frame, 
                 frame_path=frame_path, 
                 surface_beadnames=surface_beadnames, 
                 continuity_dist=continuity_dist, 
                 coarse_res=coarse_res, 
                 fine_res=fine_res)

cpus_available = os.cpu_count() - 1 if os.cpu_count() > 1 else 1
max_workers = min(cpus_available, len(frame_list))
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    for result in executor.map(worker, idx_list):
        idx, upper_area, lower_area, avg_area, upper_apl_val, lower_apl_val, avg_apl_val = result
        upper_areas[idx] = upper_area
        lower_areas[idx] = lower_area
        avg_areas[idx] = avg_area
        upper_apl[idx] = upper_apl_val
        lower_apl[idx] = lower_apl_val
        avg_apl[idx] = avg_apl_val

_ = make_outdir(outpath)
area_time_series = membpy.plot_time_series({'Upper Leaflet': upper_areas, 'Lower Leaflet': lower_areas, 'Average': avg_areas}, 
                                                 time_list, 
                                                 'Undulation Area Over Time',
                                                 'Time (ns)', 'Undulation Area (nm^2)')
area_time_series.write_html(os.path.join(outpath, 'undulation_area_time_series.html'))
apl_time_series = membpy.plot_time_series({'Upper Leaflet': upper_apl, 'Lower Leaflet': lower_apl, 'Average': avg_apl},
                                          time_list,
                                          'Area per Lipid Over Time',
                                          'Time (ns)', 'Area per Lipid (nm^2)')
apl_time_series.write_html(os.path.join(outpath, 'area_per_lipid_time_series.html'))

area_boxplot = membpy.box_plot({'Upper Leaflet': upper_areas, 'Lower Leaflet': lower_areas, 'Average': avg_areas},
                               'Distribution of Undulation Area Measurements',
                               xlabel=None, ylabel='Undulation Area (nm^2)')
area_boxplot.write_html(os.path.join(outpath, 'undulation_area_boxplot.html'))
apl_pdf = membpy.plot_pdf({'Upper Leaflet': upper_apl, 'Lower Leaflet': lower_apl, 'Average': avg_apl},
                          'Distribution of Area per Lipid Measurements',
                          xlabel=None, ylabel='Area per Lipid (nm^2)')
apl_pdf.write_html(os.path.join(outpath, 'area_per_lipid_pdf.html'))