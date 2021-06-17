"""
Some test cases to make sure we can specify various types of inputs to
the waq engine for DFM, and to test implementing new processes
"""
import os
import six
import subprocess

import stompy.model.hydro_model as hm
import stompy.model.delft.dflow_model as dfm
import stompy.model.delft.io as dio
from stompy.spatial import field

from stompy.grid import unstructured_grid
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


##
import local_config
six.moves.reload_module(local_config)

class RectangleDomain(local_config.LocalConfig,dfm.DFlowModel):
    Ly=60
    Lx=1000
    dx=20
    dy=20
    u_avg=0.25
    eta_out=0.5
    z_bed=-2.0
    
    def __init__(self,*a,**k):
        super(RectangleDomain,self).__init__(*a,**k)
        self.set_grid()
        self.mdu['output','MapFormat']=4
        self.mdu['output','MapInterval']=3600
        self.dwaq=dfm.WaqModel(self)
        self.set_bcs()

    def set_bcs(self):
        self.set_inflow()
        self.set_stage()
        
    def set_inflow(self):
        Q=self.u_avg * (self.Ly*(self.eta_out - self.z_bed))
        
        bc=hm.FlowBC(flow=Q,geom=np.array([ [0,0], [0,self.Ly]]),
                     name="inflow")
        self.add_bcs([bc])
        return bc

    def set_stage(self):
        bc=hm.StageBC(water_level=self.eta_out,
                      geom=np.array([ [self.Lx,0],[self.Lx,self.Ly]]),
                      name="stage")
        self.add_bcs([bc])
        return bc
        
    def write(self):
        # Do this first, as it updates
        # a few things in self.mdu
        self.dwaq.write_waq()
        
        super(RectangleDomain,self).write()
        
    def set_grid(self):
        g=unstructured_grid.UnstructuredGrid(max_sides=4)
        g.add_rectilinear(p0=[0,0],
                          p1=[self.Lx,self.Ly],
                          nx=1+self.Lx//self.dx,
                          ny=1+self.Ly//self.dy)
        
        g.add_node_field('node_z_bed',self.z_bed * np.ones(g.Nnodes()))
        super(RectangleDomain,self).set_grid(g)

    def run_simulation(self,threads=1,extra_args=[]):
        """
        Add in the Delwaq command line arguments
        """
        extra=extra_args+["--processlibrary",self.waq_proc_def]
        return super(RectangleDomain,self).run_simulation(threads=threads,
                                                          extra_args=extra)

class AgeDomain(RectangleDomain):
    def __init__(self,*a,**k):
        super(AgeDomain,self).__init__(*a,**k)

        # Basic age:
        self.dwaq.add_substance(name='NO3',active=True)
        self.dwaq.add_substance(name='RcNit',active=True)
        self.dwaq.add_param(name="TcNit",value=1)
        self.dwaq.add_param(name="NH4",value=1)
        self.dwaq.add_process(name="Nitrif_NH4")
        
    def set_inflow(self):
        flow_bc=super(AgeDomain,self).set_inflow()
        conc_bc=dfm.DelwaqScalarBC(parent=flow_bc,scalar='RcNit',value=1)
        self.add_bcs([conc_bc])

class PartialAgeDomain(RectangleDomain):
    def __init__(self,*a,**k):
        super(PartialAgeDomain,self).__init__(*a,**k)

        # Basic age:
        self.dwaq.add_substance(name='NO3',active=True)
        self.dwaq.add_substance(name='RcNit',active=True)
        self.dwaq.add_param(name="TcNit",value=1)
        self.dwaq.add_process(name="Nitrif_NH4")

        # For partial age, use NH4 to select spatially
        # This is supposed to end up in the ext file, such as waqparameterNH4
        # Could have different types for value, and if it's a type that indicates
        # spatial parameter, then it gets different handling in write_param() ?
        # Manual says that a spatial or temporal variable set via .ext file
        # overrules a constant added to the substances file.
        self.dwaq.add_param(name="NH4",value=1)

    def write_forcing(self):
        super(PartialAgeDomain,self).write_forcing()
        #self.write_partial_poly()
        #self.write_partial_raster()
        self.write_partial_constant_raster()
        
    def write_partial_poly(self):
        """
        Tests specifying NH4 by polygon in order to get 
        partial age
        """
        bc_fn=self.ext_force_file()
        with open(bc_fn,'at') as fp:
            stanza=f"""
QUANTITY=waqparameterNH4
FILENAME=partial.pol
FILETYPE=10
METHOD=4
OPERAND=O
VALUE=1
"""
            fp.write(stanza)
        pli_fn=os.path.join(self.run_dir,'partial.pol')
        dio.write_pli(pli_fn, [ ('partial',np.array([ [400,0],
                                                      [600,0],
                                                      [600,50],
                                                      [400,50] ])) ])

    def write_partial_constant_raster(self):
        bc_fn=self.ext_force_file()
        with open(bc_fn,'at') as fp:
            # filetype 11: netcdf gridded data
            # method 3: interpolate time and space, save coefficients
            # How to structure the nc file?
            # There is one example file in the distribution
            # "04_testcases/Example Testcases/f011_wind/c055_curvi_netcdf/wind_deltares_matroos_hirlam72_20150418.nc"

            stanza=f"""
QUANTITY=waqsegmentfunctionNH4
FILENAME=partial.nc
VARNAME=nh4
FILETYPE=11
METHOD=3
OPERAND=O
"""
            fp.write(stanza)

        # Use a field to define the raster
        nh4=field.SimpleGrid( extents=[0,self.Lx,0,self.Ly],
                              F=np.zeros([100,30]) )

        X,Y=nh4.XY()
        nh4.F=1.0*(X>550)
        
        x,y=nh4.xy()

        ds=xr.Dataset()
        ds['x']=('x',), x
        ds['x'].attrs['standard_name']='projection_x_coordinate'
        ds['x'].attrs['units']='m'
        ds['y']=('y',), y
        ds['y'].attrs['standard_name']='projection_y_coordinate'
        ds['y'].attrs['units']='m'

        ds['nh4']=('time','y','x'), np.stack([nh4.F,nh4.F])
        # Must be at least 2 timesteps
        ds['time']=('time',), np.r_[self.run_start,self.run_stop]

        ds.to_netcdf( os.path.join(self.run_dir,'partial.nc') )

    def set_inflow(self):
        flow_bc=super(PartialAgeDomain,self).set_inflow()
        conc_bc=dfm.DelwaqScalarBC(parent=flow_bc,scalar='RcNit',value=1)
        self.add_bcs([conc_bc])

def test_partial_age():
    # Age accumulates only in right/downstream half of the domain.
    model=PartialAgeDomain(run_start=np.datetime64("2016-01-01 00:00"),
                           run_stop =np.datetime64("2016-01-04 00:00"),
                           run_dir="run_partial")

    model.write()
    model.partition()
    try:
        model.run_simulation()
    except subprocess.CalledProcessError as exc:
        print("Run failed, with output:")
        print(exc.output)
        raise

    maps=model.map_outputs()
    ds=xr.open_dataset(maps[0])
    ds.close()
    ds=xr.open_dataset(maps[0])

    g=unstructured_grid.UnstructuredGrid.read_ugrid(ds)

    ubar=ds['mesh2d_ucmag'].isel(time=-1).values.mean()
    age=ds['mesh2d_NO3'].isel(time=-1).values / ds['mesh2d_RcNit'].isel(time=-2).values
    # loose validation
    age_max_expected=(model.Lx/2) / ubar / 86400.
    rel_err=np.abs(age.max() - age_max_expected) / age_max_expected
    assert rel_err<0.02,"Partial age doesn't look right"
    left_age=age[ g.cells_center()[:,0]<model.Lx/2 ]
    rel_err=np.abs(left_age/age_max_expected)
    assert rel_err.max()<0.01,"Shouldn't see age on the left"

    fig=plt.figure(1)
    fig.clf()

    scalars=['mesh2d_ucmag','mesh2d_RcNit','age']

    fig,axs=plt.subplots(len(scalars),1,sharex=True,num=1)

    for ax,scal in zip(axs,scalars):
        g.plot_edges(color='k',lw=0.5,ax=ax,)

        if scal=='age':
            values=ds['mesh2d_NO3'].isel(time=-1).values / ds['mesh2d_RcNit'].isel(time=-2).values
        else:
            values=ds[scal].isel(time=-1)
        ccoll=g.plot_cells(values=values,cmap='jet',ax=ax)

        ax.axis('equal')
        plt.colorbar(ccoll,ax=ax,label=scal)
    fig.savefig(os.path.join(model.run_dir,"age-plot.png"))

# test_partial_age()    
## 

class SpatioTemporalAgeDomain(PartialAgeDomain):
    def write_forcing(self):
        super(SpatioTemporalAgeDomain,self).write_forcing()
        self.write_spatiotemporal_raster()
    
    def write_spatiotemporal_raster(self):
        bc_fn=self.ext_force_file()
        with open(bc_fn,'at') as fp:
            # filetype 11: netcdf gridded data
            # method 3: interpolate time and space, save coefficients
            # How to structure the nc file?
            # There is one example file in the distribution
            # "04_testcases/Example Testcases/f011_wind/c055_curvi_netcdf/wind_deltares_matroos_hirlam72_20150418.nc"

            stanza=f"""
QUANTITY=waqsegmentfunctionNH4
FILENAME=spatiotemporal.nc
VARNAME=nh4
FILETYPE=11
METHOD=3
OPERAND=O
"""
            fp.write(stanza)

        # Use a field to define the raster
        nh4=field.SimpleGrid( extents=[0,self.Lx,0,self.Ly],
                              F=np.zeros([100,30]) )

        X,Y=nh4.XY()
        nh4.F=1.0*(X>550)
        
        x,y=nh4.xy()

        ds=xr.Dataset()
        ds['x']=('x',), x
        ds['x'].attrs['standard_name']='projection_x_coordinate'
        ds['x'].attrs['units']='m'
        ds['y']=('y',), y
        ds['y'].attrs['standard_name']='projection_y_coordinate'
        ds['y'].attrs['units']='m'

        ds['nh4']=('time','y','x'), np.stack([0*nh4.F,0*nh4.F,nh4.F])
        t_mid = self.run_start + (self.run_stop - self.run_start)/2
        ds['time']=('time',), np.r_[self.run_start,t_mid,self.run_stop]

        ds.to_netcdf( os.path.join(self.run_dir,'spatiotemporal.nc') )

#def test_spatiotemporal_age():
if 1:
    # Age accumulates only in right/downstream half of the domain,
    # and only during the second half of the run
    model=SpatioTemporalAgeDomain(run_start=np.datetime64("2016-01-01 00:00"),
                                  run_stop =np.datetime64("2016-01-04 00:00"),
                                  run_dir="run_spatiotemporal")

    model.write()
    model.partition()
    try:
        model.run_simulation()
    except subprocess.CalledProcessError as exc:
        print("Run failed, with output:")
        print(exc.output)
        raise

    maps=model.map_outputs()
    ds=xr.open_dataset(maps[0])
    ds.close()
    ds=xr.open_dataset(maps[0])

    
        
# Basic age:
# Specify a constant in time, spatially variable.

# At 0.25 m/s, domain is 1000m long, so max age should be (1000/0.25) / 86400
# or 0.046.  That checks out.

##

fig=plt.figure(1)

scalars=['mesh2d_ucmag','mesh2d_RcNit','age']


for tidx in range(0,ds.dims['time'],10):
    print(tidx)
    fig.clf()
    fig,axs=plt.subplots(len(scalars),1,sharex=True,num=1)
    
    for ax,scal in zip(axs,scalars):
        g.plot_edges(color='k',lw=0.5,ax=ax,)

        if scal=='age':
            values=ds['mesh2d_NO3'].isel(time=tidx).values / ds['mesh2d_RcNit'].isel(time=tidx).values
        else:
            values=ds[scal].isel(time=tidx)
        ccoll=g.plot_cells(values=values,cmap='jet',ax=ax)

        ax.axis('equal')
        plt.colorbar(ccoll,ax=ax,label=scal)
    fig.canvas.draw()
    plt.pause(0.1)

# This code is working
# Can specify a polygon-based age field, raster constant, or raster/changing.

