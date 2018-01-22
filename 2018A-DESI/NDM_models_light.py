from NDM_utils import *

class DESI_NDM(object):
    def __init__(self):
        """
        Import the intersection training data. Set global parameters to default values.
        """
        # ---- Warning: These parameters are frozen and should not be changed ---- #
        # Grid parameters
        self.var_x_limits = [-.25, 3.5] # g-z
        self.var_y_limits = [-1, 1.4] # g-r
        self.gmag_limits = [19.5, 24.]
        self.num_bins = [375, 240, 450]

        # Cell_number in selection. Together with grid parameters this
        # is a representation of the selection region.
        self.cell_select = None

        return

    def var_reparam(self, gflux, rflux, zflux, oii = None):
        """
        Given the input variables return the model3 parametrization as noted above.
        """
        mu_g = flux2asinh_mag(gflux, band = "g")
        mu_r = flux2asinh_mag(rflux, band = "r")
        mu_z = flux2asinh_mag(zflux, band = "z")
        if oii is not None:
            mu_oii = flux2asinh_mag(oii, band = "oii")
            return mu_g - mu_z, mu_g - mu_r, mu_g - mu_oii, flux2mag(gflux)
        else:
            return mu_g - mu_z, mu_g - mu_r, None, flux2mag(gflux)


    def apply_selection(self, gflux, rflux, zflux):
        """
        Model 3
        Given gflux, rflux, zflux of samples, return a boolean vector that gives the selection.
        """
        mu_g = flux2asinh_mag(gflux, band = "g")
        mu_r = flux2asinh_mag(rflux, band = "r")
        mu_z = flux2asinh_mag(zflux, band = "z")

        var_x = mu_g - mu_z
        var_y = mu_g - mu_r
        gmag = flux2mag(gflux)

        samples = [var_x, var_y, gmag]

        # Generate cell number 
        cell_number = multdim_grid_cell_number(samples, 3, [self.var_x_limits, self.var_y_limits, self.gmag_limits], self.num_bins)

        # Sort the cell number
        idx_sort = cell_number.argsort()
        cell_number = cell_number[idx_sort]

        # Placeholder for selection vector
        iselect = check_in_arr2(cell_number, self.cell_select)

        # The last step is necessary in order for iselect to have the same order as the input sample variables.
        idx_undo_sort = idx_sort.argsort()        
        return iselect[idx_undo_sort]

    def cell_select_centers(self):
        """
        Return selected cells centers given cell numbers.
        """
        limits = [self.var_x_limits, self.var_y_limits, self.gmag_limits]
        Ncell_select = self.cell_select.size # Number of cells in the selection
        centers = [None, None, None]

        for i in range(3):
            Xmin, Xmax = limits[i]
            bin_edges, dX = np.linspace(Xmin, Xmax, self.num_bins[i]+1, endpoint=True, retstep=True)
            bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.
            idx = (self.cell_select % np.multiply.reduce(self.num_bins[i:])) //  np.multiply.reduce(self.num_bins[i+1:])
            centers[i] = bin_centers[idx.astype(int)]

        return np.asarray(centers).T

    def gen_select_boundary_slices(self, slice_dir = 2, save_dir="../figures/", \
        prefix = "test", output_sparse=True, increment=10, centers=None, plot_ext=False,\
        gflux_ext=None, rflux_ext=None, zflux_ext=None, ibool_ext = None,\
        var_x_ext=None, var_y_ext=None, gmag_ext=None, use_parameterized_ext=False,\
        pt_size=10, pt_size_ext=10, alpha_ext=0.5, guide=False):
        """
        Given slice direction, generate slices of boundary

        0: var_x
        1: var_y
        2: gmag

        If plot_ext True, then plot user supplied external objects.

        If centers is not None, then use it instead of generating one.

        If use_parameterized_ext, then the user may provide already parameterized version of the external data points.

        If guide True, then plot the guide line.

        If output_sparse=True, then only 10% of the boundaries are plotted and saved.
        """

        slice_var_tag = ["mu_gz", "mu_gr", "gmag"]
        var_names = ["$\mu_g - \mu_z$", "$\mu_g - \mu_r$", "$g$"]

        if centers is None:
            centers = self.cell_select_centers()

        if guide:
            x_guide, y_guide = gen_guide_line()

        if plot_ext:
            if use_parameterized_ext:
                if ibool_ext is not None:
                    var_x_ext = var_x_ext[ibool_ext]
                    var_y_ext = var_y_ext[ibool_ext]
                    gmag_ext = gmag_ext[ibool_ext]                
            else:     
                if ibool_ext is not None:
                    gflux_ext = gflux_ext[ibool_ext]
                    rflux_ext = rflux_ext[ibool_ext]
                    zflux_ext = zflux_ext[ibool_ext]

                mu_g, mu_r, mu_z = flux2asinh_mag(gflux_ext, band="g"), flux2asinh_mag(rflux_ext, band="r"), flux2asinh_mag(zflux_ext, band="z")
                var_x_ext = mu_g-mu_z
                var_y_ext = mu_g-mu_r
                gmag_ext = flux2mag(gflux_ext)

            variables = [var_x_ext, var_y_ext, gmag_ext]

        limits = [self.var_x_limits, self.var_y_limits, self.gmag_limits]        
        Xmin, Xmax = limits[slice_dir]
        bin_edges, dX = np.linspace(Xmin, Xmax, self.num_bins[slice_dir]+1, endpoint=True, retstep=True)

        print slice_var_tag[slice_dir]
        if output_sparse:
            iterator = range(0, self.num_bins[slice_dir], increment)
        else:
            iterator = range(self.num_bins[slice_dir])

        for i in iterator: 
            ibool = (centers[:, slice_dir] < bin_edges[i+1]) & (centers[:, slice_dir] > bin_edges[i])
            centers_slice = centers[ibool, :]
            fig = plt.figure(figsize=(7, 7))
            idx = range(3)
            idx.remove(slice_dir)
            plt.scatter(centers_slice[:,idx[0]], centers_slice[:,idx[1]], edgecolors="none", c="green", alpha=0.5, s=pt_size)
            if plot_ext:
                ibool = (variables[slice_dir] < bin_edges[i+1]) & (variables[slice_dir] > bin_edges[i])
                plt.scatter(variables[idx[0]][ibool], variables[idx[1]][ibool], edgecolors="none", c="red", s=pt_size_ext, alpha=alpha_ext)
            plt.xlabel(var_names[idx[0]], fontsize=25)
            plt.ylabel(var_names[idx[1]], fontsize=25)

            if guide and (slice_dir==2):
                plt.plot(x_guide, y_guide, c="orange", lw = 2)
            # plt.axis("equal")
            plt.xlim(limits[idx[0]])
            plt.ylim(limits[idx[1]])
            title_str = "%s [%.3f, %.3f]" % (var_names[slice_dir], bin_edges[i], bin_edges[i+1])
            print i, title_str
            plt.title(title_str, fontsize=25, y =1.05)
            plt.savefig(save_dir+prefix+"-boundary-%s-%d.png" % (slice_var_tag[slice_dir], i), bbox_inches="tight", dpi=200)
            plt.close()        


