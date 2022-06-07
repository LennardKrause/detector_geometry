import numpy as np
from matplotlib.widgets import RadioButtons, Slider
from matplotlib import pyplot as plt
from matplotlib import cm, colors, patches
    
def get_specs():
    ######################
    # Setup the geometry #
    ######################
    geo = container()
    geo.det_type = 'Eiger2 CdTe'   # [str]  Pilatus3 / Eiger2
    geo.det_size = '4M'            # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
    geo.ener = 21.0                # [keV]  Beam energy
    geo.dist = 75.0                # [mm]   Detector distance
    geo.yoff = 0.0                 # [mm]   Detector offset (vertical)
    geo.rota = 25.0                # [deg]  Detector rotation
    geo.tilt = 0.0                 # [deg]  Detector tilt
    geo.unit = 1                   # [0-3]  Contour legend
                                   #          0: 2-Theta
                                   #          1: d-spacing
                                   #          2: q-space
                                   #          3: sin(theta)/lambda
    geo.std_idx = 0                # [0-3]  Plot standard contours
                                   #          0: None
                                   #          1: LaB6
                                   #          2: Si
                                   #          3: CeO2
    # What standards should be available
    # The d spacings will be imported from pyFAI
    # and we use 2 lists whose ORDER must match!
    #  - this is the display name
    geo.std_names = ['None', r'$LaB_6$', r'$Si$', r'$CeO_2$']
    #  - this is what pyFAI understands
    geo.std_pyFAI = ['None', 'LaB6', 'Si', 'CeO2']

    ###########################
    # Detector Specifications #
    ###########################
    det = container()
    if geo.det_type.startswith('Pilatus'):
        ###############################
        # Specifications for Pilatus3 #
        ###############################
        det.hms = 83.8    # [mm]  Module size (horizontal)
        det.vms = 33.5    # [mm]  Module size (vertical)
        det.pxs = 172e-3  # [mm]  Pixel size
        det.hgp = 7       # [pix] Gap between modules (horizontal)
        det.vgp = 17      # [pix] Gap between modules (vertical)
        det.cbh = 0       # [mm]  Central beam hole
        det.name = f'{geo.det_type} {geo.det_size}'
        det.sizes = {'300K':(1,3),'1M':(2,5),'2M':(3,8),'6M':(5,12)}
        if geo.det_size not in det.sizes.keys():
            print('Unknown detector type/size combination!')
            raise SystemExit
        det.hmn, det.vmn = det.sizes[geo.det_size]
    elif geo.det_type.startswith('Eiger'):
        #############################
        # Specifications for Eiger2 #
        #############################
        det.hms = 77.1    # [mm]  Module size (horizontal)
        det.vms = 38.4    # [mm]  Module size (vertical)
        det.pxs = 75e-3   # [mm]  Pixel size
        det.hgp = 38      # [pix] Gap between modules (horizontal)
        det.vgp = 12      # [pix] Gap between modules (vertical)
        det.cbh = 0       # [mm]  Central beam hole
        det.name = f'{geo.det_type} {geo.det_size}'
        det.sizes = {'1M':(1,2),'4M':(2,4),'9M':(3,6),'16M':(4,8)}
        if geo.det_size not in det.sizes.keys():
            print('Unknown detector type/size combination!')
            raise SystemExit
        det.hmn, det.vmn = det.sizes[geo.det_size]
    elif geo.det_type.startswith('MPCCD'):
        #############################
        # Specifications for MPCCD #
        #############################
        det.hms = 51.2    # [mm]  Module size (horizontal)
        det.vms = 25.6    # [mm]  Module size (vertical)
        det.pxs = 50e-3   # [mm]  Pixel size
        det.hgp = 18      # [pix] Gap between modules (horizontal)
        det.vgp = 27      # [pix] Gap between modules (vertical)
        det.cbh = 3       # [mm]  Central beam hole
        det.name = f'{geo.det_type} {geo.det_size}'
        det.sizes = {'4M':(2,4)}
        if geo.det_size not in det.sizes.keys():
            print('Unknown detector type/size combination!')
            raise SystemExit
        det.hmn, det.vmn = det.sizes[geo.det_size]
    else:
        ###########################################
        # ADD CUSTOM DETECTOR SPECIFICATIONS HERE #
        ###########################################
        det.hms = 100.0   # [mm]  Module size (horizontal)
        det.vms = 140.0   # [mm]  Module size (vertical)
        det.pxs = 10e-3   # [mm]  Pixel size
        det.hgp = 0       # [pix] Gap between modules (horizontal)
        det.vgp = 0       # [pix] Gap between modules (vertical)
        det.hmn = 1       # [int] Number of modules (horizontal)
        det.vmn = 1       # [int] Number of modules (vertical)
        det.cbh = 0       # [mm]  Central beam hole
        det.name = f'{geo.det_type} {geo.det_size} Octal'

    ################
    # Plot Details #
    ################
    plo = container()
    # - geometry contour section - 
    plo.cont_tth_min = 5                # [int]    minimum 2-theta contour line
    plo.cont_tth_max = 120              # [int]    maximum 2-theta contour line
    plo.cont_tth_num = 24               # [int]    number of contour lines
    plo.cont_geom_cmark = 'o'           # [marker] Beam center marker (geometry)
    plo.cont_geom_csize = 4             # [int]    Beam center size (geometry)
    plo.cont_geom_alpha = 1.00          # [float]  Contour alpha (geometry)
    plo.cont_geom_cmap_name = 'viridis' # [cmap]   Contour colormap (geometry)
    # - standard contour section - 
    plo.cont_standard = True            # [bool]   Plot additional contour lines
                                        #          e.g. a LaB6 standard
    plo.cont_std_alpha = 0.25           # [float]  Standard contour alpha
    plo.cont_std_color = 'gray'         # [color]  Standard contour color
    plo.cont_std_lw = 2.5               # [float]  Standard contour linewidth
    # - normal incidence contour section - 
    plo.cont_norm_inc = False           # [bool]   Plot additional contour lines
                                        #          for normal incidence geometry
    plo.cont_norm_cmark = 'o'           # [marker] Beam center marker (original)
    plo.cont_norm_csize = 4             # [int]    Beam center size (original)
    plo.cont_norm_alpha = 0.25          # [float]  Contour alpha (original)
    plo.cont_norm_color = 'gray'        # [color]  Contour color (original)
    # - general section - 
    plo.cont_reso_min = 50              # [int]    Minimum contour steps
    plo.cont_reso_max = 500             # [int]    Maximum contour steps
    plo.module_alpha = 0.20             # [float]  Detector module alpha
    plo.module_color = 'gray'           # [color]  Detector module color
    plo.margin_top = 0.95               # [float]  Plot margin for title
    plo.plot_size = 8                   # [int]    Plot size
    plo.label_size = 9                  # [int]    Label size
    plo.plot_dpi = 300                  # [int]    Set plot DPI for saving
    plo.plot_color = 0.35               # [float]  Button color from colormap (0.0 - 1.0)
                                        # [str]    Button color e.g. '#1f77b4'
    plo.interactive = True              # [bool]   Make the plot interactive
    plo.action_ener = True              # [bool]   Show energy slider
    plo.action_dist = True              # [bool]   Show distance slider
    plo.action_rota = True              # [bool]   Show rotation slider
    plo.action_yoff = True              # [bool]   Show offset slider
    plo.action_tilt = True              # [bool]   Show tilt slider
    plo.action_radio = True             # [bool]   Show radio buttons

    ##########
    # Limits #
    ##########
    lmt = container()
    lmt.ener_min = 1.0   # [float] Energy minimum [keV]
    lmt.ener_max = 100.0 # [float] Energy maximum [keV]
    lmt.ener_stp = 1.0   # [float] Energy step size [keV]
    lmt.dist_min = 40.0  # [float] Distance minimum [mm]
    lmt.dist_max = 150.0 # [float] Distance maximum [mm]
    lmt.dist_stp = 1.0   # [float] Distance step size [mm]
    lmt.yoff_min = 0.0   # [float] Offset minimum [mm]
    lmt.yoff_max = 200.0 # [float] Offset maximum [mm]
    lmt.yoff_stp = 1.0   # [float] Offset step size [mm]
    lmt.rota_min = 0.0   # [float] Rotation minimum [deg]
    lmt.rota_max = 75.0  # [float] Rotation maximum [deg]
    lmt.rota_stp = 1.0   # [float] Rotation step size [deg]
    lmt.tilt_min = 0.0   # [float] Tilt minimum [deg]
    lmt.tilt_max = 45.0  # [float] Tilt maximum [deg]
    lmt.tilt_stp = 1.0   # [float] Tilt step size [deg]

    ###################################
    # !!! Don't change below here !!! #
    ###################################
    return geo, det, plo, lmt

def main():
    # fetch the geometry, detector, plot specifications and limits
    geo, det, plo, lmt = get_specs()
    # translate unit for plot title
    geo.unit_names = [r'2$\theta$', r'$d_{space}$', r'$q_{space}$', r'$sin(\theta)/\lambda$']
    if geo.unit >= len(geo.unit_names):
        print(f'Error: Valid geo.unit range is from 0 to {len(geo.unit_names)-1}, geo.unit={geo.unit}')
        raise SystemExit
    # get colormap from name
    plo.cont_geom_cmap = cm.get_cmap(plo.cont_geom_cmap_name)
    # figure out the color of the buttons and slider handles
    try:
        # try to derive color from colormap
        plo.plot_handle_color = plo.cont_geom_cmap(plo.plot_color)
    except TypeError:
        # use color as defined by user
        plo.plot_handle_color = plo.plot_color
    # import pyFAI if standard contours are enabled
    if plo.cont_standard:
        from pyFAI import calibrant
        geo.pyFAI_calibrant = calibrant
        # get contour lines f contours are already selected (index is not 0, not None)
        if geo.std_idx > 0:
            # get the d spacings for the calibrtant from pyFAI
            plo.cont_std_dsp = np.array(geo.pyFAI_calibrant.get_calibrant(geo.std_pyFAI[geo.std_idx]).get_dSpacing())
    # set rcParams
    plt.rcParams['savefig.dpi'] = plo.plot_dpi
    # figure out proper plot dimensions
    plo.xdim = (det.hms * det.hmn + det.pxs * det.hgp * det.hmn + det.cbh)/2
    plo.ydim = (det.vms * det.vmn + det.pxs * det.vgp * det.vmn + det.cbh)/2
    plo.fig_ratio = plo.xdim / plo.ydim
    # scale contour grid to detector size
    plo.cont_grid_max = int(np.ceil(max(plo.xdim, plo.ydim)))
    # generate contour levels
    plo.cont_levels = np.linspace(plo.cont_tth_min, plo.cont_tth_max, plo.cont_tth_num)
    # init the plot
    fig = plt.figure()
    # add axes for the detector modules
    bg = fig.add_subplot(111, aspect='equal')
    # add axes for the contours
    ax = fig.add_subplot(111, aspect='equal')
    # remove both axis, ticks and labels
    bg.set_axis_off()
    ax.set_axis_off()
    # limit the axis x and y
    bg.set_xlim(-plo.xdim, plo.xdim)
    ax.set_xlim(-plo.xdim, plo.xdim)
    bg.set_ylim(-plo.ydim, plo.ydim)
    ax.set_ylim(-plo.ydim, plo.ydim)
    # setup detector and geometry
    build_detector(bg, det, plo)
    # create cones and draw contour lines
    draw_contours(ax, geo, plo)
    # generate some sense of interactivity
    if plo.interactive:
        # add short title to make room for the text boxes
        plt.suptitle(f'{det.name}', size=10, fontweight='bold')
        # define sliders
        if plo.action_ener:
            # make room for the sliders
            plo.margin_top -= 0.02
            plo.margin_right = 0.6
            # make space for the calibrant picker
            if plo.cont_standard:
                plo.margin_right -= 0.1
            # add slider
            sli_ener = add_slider('Energy [keV] ' , 'ener', 0.3, plo.margin_top, plo.margin_right, 0.025, geo.ener, lmt.ener_min, lmt.ener_max, lmt.ener_stp, fig, ax, geo, plo)
        if plo.action_dist:
            plo.margin_top -= 0.02
            sli_dist = add_slider('Distance [mm] ', 'dist', 0.3, plo.margin_top, plo.margin_right, 0.025, geo.dist, lmt.dist_min, lmt.dist_max, lmt.dist_stp, fig, ax, geo, plo)
        if plo.action_yoff:
            plo.margin_top -= 0.02
            sli_yoff = add_slider('Offset [mm] '  , 'yoff', 0.3, plo.margin_top, plo.margin_right, 0.025, geo.yoff, lmt.yoff_min, lmt.yoff_max, lmt.yoff_stp, fig, ax, geo, plo)
        if plo.action_tilt:
            plo.margin_top -= 0.02
            sli_tilt = add_slider('Tilt [˚] '     , 'tilt', 0.3, plo.margin_top, plo.margin_right, 0.025, geo.tilt, lmt.tilt_min, lmt.tilt_max, lmt.tilt_stp, fig, ax, geo, plo)
        if plo.action_rota:
            plo.margin_top -= 0.02
            sli_rota = add_slider('Rotation [˚] ' , 'rota', 0.3, plo.margin_top, plo.margin_right, 0.025, geo.rota, lmt.rota_min, lmt.rota_max, lmt.rota_stp, fig, ax, geo, plo)
        # add radio buttons and an axis for the buttons
        if plo.action_radio:
            # figure out a proper size of the axis
            _ds = 1.0 - (plo.margin_top-0.01)
            # add the axes
            axs_unit = fig.add_axes([0.0, plo.margin_top-0.01, _ds, _ds], frameon=False, aspect='equal')
            box_unit = RadioButtons(axs_unit, geo.unit_names, active=geo.unit, activecolor=plo.plot_handle_color)
            box_unit.on_clicked(lambda val: update_plot('unit', geo.unit_names.index(val), fig, geo, plo, ax))
            # change label size
            for l in box_unit.labels:
                l.set_size(plo.label_size)
        # add radio buttons and an axis for the buttons
        if plo.cont_standard:
            # figure out a proper size of the axis
            _ds = 1.0 - (plo.margin_top-0.01)
            # add the axes
            axs_std = fig.add_axes([0.85, plo.margin_top-0.01, _ds, _ds], frameon=False, aspect='equal')
            box_std = RadioButtons(axs_std, geo.std_names, active=geo.std_idx, activecolor=plo.plot_handle_color)
            box_std.on_clicked(lambda val: update_plot('std', geo.std_names.index(val), fig, geo, plo, ax))
            # change label size
            for l in box_std.labels:
                l.set_size(plo.label_size)
    else:
        # make room for the second title line
        plo.margin_top -= 0.02
        # add title / information
        plt.suptitle(f'{det.name} | Energy: {geo.ener} keV | Distance: {geo.dist} cm\nRotation: {geo.rota}° | Tilt: {geo.tilt}° | Offset: {geo.yoff} cm | Units: {geo.unit_names[geo.unit]}', size=10)
    # adjust the margins
    plt.subplots_adjust(top=plo.margin_top, bottom=0, right=1, left=0, hspace=0, wspace=0)
    # adjust the figure size
    fig.set_size_inches(plo.plot_size * plo.margin_top * plo.fig_ratio, plo.plot_size)
    # show the plot
    plt.show()

def build_detector(bg, det, plo):
    # build detector modules
    # beam position is between the modules (even) or at the center module (odd)
    # determined by the "+det.hmn%2" part
    for i in range(-det.hmn//2+det.hmn%2, det.hmn-det.hmn//2):
        for j in range(-det.vmn//2+det.vmn%2, det.vmn-det.vmn//2):
            # - place modules along x (i) and y (j) keeping the gaps in mind ( + (det.hgp*det.pxs)/2)
            # - the " - ((det.hms+det.hgp*det.pxs)/2)" positions the origin (the beam) at the center of a module
            #   and "det.hmn%2" makes sure this is only active for detectors with an odd number of modules
            # - define sets of panels that collectively move to realize a central hole offset for MPCCD detectors
            #   that are used at SACLA/SPring-8:
            #   x = (...) + (det.cbh/2)*(2*(j&det.vmn)//det.vmn-1)
            #   y = (...) + (det.cbh/2)*(1-2*(i&det.hmn)//det.hmn)
            # - negative values of det.cbh for 'clockwise' offset order
            origin_x = i*(det.hms+det.hgp*det.pxs) - ((det.hms+det.hgp*det.pxs)/2)*(det.hmn%2) + (det.hgp*det.pxs)/2 + (det.cbh/2)*(2*(j&det.vmn)//det.vmn-1)
            origin_y = j*(det.vms+det.vgp*det.pxs) - ((det.vms+det.vgp*det.pxs)/2)*(det.vmn%2) + (det.vgp*det.pxs)/2 + (det.cbh/2)*(1-2*(i&det.hmn)//det.hmn)
            # DEBUG: print indices on panels
            #bg.annotate(f'{i} {j}', (origin_x+det.hms/2, origin_y+det.vms/2), color='gray', alpha=0.2, size=16, ha='center')
            # add the module
            bg.add_patch(patches.Rectangle((origin_x, origin_y),  det.hms, det.vms, color=plo.module_color, alpha=plo.module_alpha))

def draw_contours(ax, geo, plo):
    # calculate the offset of the contours resulting from yoff and rotation
    # shift the grid to draw the cones, to make sure the contours are drawn
    # within the visible area
    _comp_shift = -(geo.yoff + np.tan(np.deg2rad(geo.rota))*geo.dist)
    # increase the the cone grid to allow more
    # contours to be drawn as the plane is tilted
    _comp_add = np.tan(np.deg2rad(geo.tilt))*geo.dist
    # draw beam center
    if plo.cont_norm_inc:
        ax.plot(0, 0, color=plo.cont_norm_color, marker=plo.cont_norm_cmark, ms=plo.cont_norm_csize, alpha=plo.cont_norm_alpha)
    ax.plot(0, _comp_shift, color=colors.to_hex(plo.cont_geom_cmap(1)), marker=plo.cont_geom_cmark, ms=plo.cont_geom_csize, alpha=plo.cont_geom_alpha)
    # draw contour lines
    for _n,_ttd in enumerate(plo.cont_levels):
        # convert theta in degrees to radians
        _ttr = np.deg2rad(_ttd)
        # calculate ratio of sample to detector distance (sdd)
        # and contour distance to beam center (cbc)
        # _rat = sdd/cbc = 1/tan(2-theta)
        # this is used to scale the cones Z dimension
        _rat = 1/np.tan(_ttr)
        # apply the min/max grid resolution
        _grd_res = max(min(int(plo.cont_reso_min*_rat),plo.cont_reso_max), plo.cont_reso_min)
        # prepare the grid for the cones/contours
        # adjust the resolution using i (-> plo.cont_levels),
        # as smaller cones/contours (large i) need higher sampling
        # but make sure the sampling rate doesn't fall below the
        # user set plo.cont_reso_min value and plo.cont_reso_max
        # prevents large numbers that will take seconds to draw
        _x0 = np.linspace(-plo.cont_grid_max, plo.cont_grid_max, _grd_res)
        # the grid position needs to adjusted upon change of geometry
        # the center needs to be shifted by _geo_offset to make sure 
        # sll contour lines are drawn
        _x1 = np.linspace(-plo.cont_grid_max + _comp_shift, plo.cont_grid_max - _comp_shift + _comp_add, _grd_res)
        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        _stl = np.sin(_ttr/2)/(12.398/geo.ener)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        _dsp = 1/(2*_stl)
        # prepare the values in the different units / labels
        _units = {0:np.rad2deg(_ttr), 1:_dsp, 2:_stl*4*np.pi, 3:_stl}
        # draw additional contours for normal incidence geometry
        if plo.cont_norm_inc:
            X0, Y0 = np.meshgrid(_x0,_x0)
            Z0 = np.sqrt(X0**2+Y0**2)*_rat
            X,Y,Z = geo_cone(X0, Y0, Z0, 0, 0, 0, geo.dist)
            # don't draw contour lines that are out of bounds
            # make sure Z is large enough to draw the contour
            if np.max(Z) >= geo.dist:
                c0 = ax.contour(X, Y, Z, [geo.dist], colors=plo.cont_norm_color, alpha=plo.cont_norm_alpha)
                # label original geometry contours
                fmt = {c0.levels[0]:f'{np.round(_units[geo.unit],2):.2f}'}
                ax.clabel(c0, c0.levels, inline=True, fontsize=plo.label_size, fmt=fmt, manual=[(plo.xdim,plo.ydim)])
        # draw contours for the tilted/rotated/moved geometry
        # use the offset adjusted value x1 to prepare the grid
        X0, Y0 = np.meshgrid(_x1,_x0)
        Z0 = np.sqrt(X0**2+Y0**2)*_rat
        X,Y,Z = geo_cone(X0, Y0, Z0, geo.rota, geo.tilt, geo.yoff, geo.dist)
        # make sure Z is large enough to draw the contour
        if np.max(Z) > geo.dist:
            c1 = ax.contour(X, Y, Z, [geo.dist], colors=colors.to_hex(plo.cont_geom_cmap((_n+1)/len(plo.cont_levels))), alpha=plo.cont_geom_alpha)
            # label moved geometry contours
            fmt = {c1.levels[0]:f'{np.round(_units[geo.unit],2):.2f}'}
            ax.clabel(c1, c1.levels, inline=True, fontsize=plo.label_size, fmt=fmt, manual=[(0,plo.ydim)])
        else:
            # if the Z*i is too small, break as the following cycles
            # will make Z only smaller -> leaving no contours to draw
            # - only True if the contours are iterated high to low!
            break
    # plot standard contour lines
    if plo.cont_standard and geo.std_idx > 0:
        # this assumes that the last cycle runs for the highest resolution
        # so _dsp holds the correct maximum value up to which the 
        # satndard contour lines are to be drawn
        for _d in plo.cont_std_dsp[plo.cont_std_dsp > _dsp]:
            # lambda = 2 * d * sin(theta)
            # 2-theta = 2 * (lambda / 2*d)
            # lambda -> (12.398/geo_energy)
            _ttr = 2 * np.arcsin((12.398/geo.ener) / (2*_d))
            # calculate ratio of sample to detector distance (sdd)
            # and contour distance to beam center (cbc)
            # _rat = sdd/cbc = 1/tan(2-theta)
            # this is used to scale the cones Z dimension
            _rat = 1/np.tan(_ttr)
            # apply the min/max grid resolution
            _grd_res = max(min(int(plo.cont_reso_min*_rat),plo.cont_reso_max), plo.cont_reso_min)
            # prepare the grid for the cones/contours
            # adjust the resolution using i (-> plo.cont_levels),
            # as smaller cones/contours (large i) need higher sampling
            # but make sure the sampling rate doesn't fall below the
            # user set plo.cont_reso_min value and plo.cont_reso_max
            # prevents large numbers that will take seconds to draw
            _x0 = np.linspace(-plo.cont_grid_max, plo.cont_grid_max, _grd_res)
            # the grid position needs to adjusted upon change of geometry
            # the center needs to be shifted by _geo_offset to make sure 
            # sll contour lines are drawn
            _x1 = np.linspace(-plo.cont_grid_max + _comp_shift, plo.cont_grid_max - _comp_shift + _comp_add, _grd_res)
            # draw contours for the tilted/rotated/moved geometry
            # use the offset adjusted value x1 to prepare the grid
            X0, Y0 = np.meshgrid(_x1,_x0)
            Z0 = np.sqrt(X0**2+Y0**2)*_rat
            X,Y,Z = geo_cone(X0, Y0, Z0, geo.rota, geo.tilt, geo.yoff, geo.dist)
            # make sure Z is large enough to draw the contour
            if np.max(Z) > geo.dist:
                c1 = ax.contour(X, Y, Z, [geo.dist], colors=plo.cont_std_color, alpha=plo.cont_std_alpha, linewidths=plo.cont_std_lw)

def geo_cone(X, Y, Z, rota, tilt, yoff, dist):
    # combined rotation, tilt 'movement' is compensated
    a = np.deg2rad(tilt) + np.deg2rad(rota)
    # rotate the sample around y
    t = np.transpose(np.array([X,Y,Z]), (1,2,0))
    # rotation matrix
    m = [[np.cos(a), 0, np.sin(a)],[0,1,0],[-np.sin(a), 0, np.cos(a)]]
    # apply rotation
    X,Y,Z = np.transpose(np.dot(t, m), (2,0,1))
    # compensate for tilt not rotating
    # - revert the travel distance
    comp = np.deg2rad(tilt) * dist
    return Y,X+comp-yoff,Z

def add_slider(label, name, left, bottom, width, height, val, vmin, vmax, step, fig, ax, geo, plo):
    axs = fig.add_axes([left, bottom, width, height])
    sli = Slider(axs, label, valmin=vmin, valmax=vmax, valinit=val, handle_style={'size':plo.label_size}, valstep=step, color=plo.plot_handle_color)
    sli.vline.set_alpha(0) # Remove the mark on the slider
    sli.on_changed(lambda val: update_plot(name, val, fig, geo, plo, ax))
    sli.label.set_size(plo.label_size)
    sli.valtext.set_size(plo.label_size)
    return sli

def update_plot(nam, val, fig, geo, plo, ax):
    ##################################################
    # This is a sloppy and hacky way to achieve some #
    #   interactivity without building a proper GUI  #
    ##################################################
    # this works with MacOS backend
    # but seems to have problems with
    # Tk and Qt Agg backends for which the
    # radio buttons have to be pressed twice
    if nam == 'dist':
        geo.dist = float(val)
    elif nam == 'rota':
        geo.rota = float(val)
    elif nam == 'tilt':
        geo.tilt = float(val)
    elif nam == 'yoff':
        geo.yoff = float(val)
    elif nam == 'unit':
        geo.unit = int(val)
    elif nam == 'ener':
        geo.ener = float(val)
    elif nam == 'std':
        geo.std_idx = int(val)
        if geo.std_idx > 0:
            # get the d spacings for the calibrtant from pyFAI
            plo.cont_std_dsp = np.array(geo.pyFAI_calibrant.get_calibrant(geo.std_pyFAI[geo.std_idx]).get_dSpacing())
    # clear the axis
    # the beam center is a line
    for _a in ax.lines:
        _a.remove()
    # the contour label is a text
    for _t in ax.texts:
        _t.remove()
    # the contour is a collection
    # and needs some special attention
    # to make sure all contours are removed
    while len(ax.collections) > 0:
        ax.collections.pop()
    # re-calculate cones and re-draw contours
    draw_contours(ax, geo, plo)
    # blit the ax
    if plt.rcParams['backend'] == 'MacOSX':
        fig.canvas.blit(ax.bbox)
    else:
        fig.canvas.draw()

class container(object):
    pass

if __name__ == '__main__':
    main()
