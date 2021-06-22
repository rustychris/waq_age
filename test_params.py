"""
Some test cases to make sure we can specify various types of inputs to
the waq engine for DFM, and to test implementing new processes
"""
import os
import six
import subprocess

import stompy.model.hydro_model as hm
import stompy.model.delft.dflow_model as dfm
import stompy.model.delft.waq_scenario as dwaq

import stompy.model.delft.io as dio
from stompy.spatial import field
from stompy import utils

from stompy.grid import unstructured_grid
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


##
import local_config
six.moves.reload_module(dwaq)
six.moves.reload_module(dfm)
six.moves.reload_module(local_config)

class RectangleDomain(local_config.LocalConfig,dfm.DFlowModel):
    Ly=60
    Lx=1000
    dx=20
    dy=20
    u_avg=0.25
    eta_out=0.5
    z_bed=-2.0
    run_start=np.datetime64("2016-01-01 00:00")
    run_stop =np.datetime64("2016-01-04 00:00")
    fig_num=1
    dwaq=True
    
    def __init__(self,*a,**k):
        super(RectangleDomain,self).__init__(*a,**k)
        self.clean_run_dir()
        self.set_grid()
        self.mdu['output','MapFormat']=4
        self.mdu['output','MapInterval']=3600
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
        
    def set_grid(self):
        g=unstructured_grid.UnstructuredGrid(max_sides=4)
        g.add_rectilinear(p0=[0,0],
                          p1=[self.Lx,self.Ly],
                          nx=1+self.Lx//self.dx,
                          ny=1+self.Ly//self.dy)
        
        g.add_node_field('node_z_bed',self.z_bed * np.ones(g.Nnodes()))
        super(RectangleDomain,self).set_grid(g)

    def run_all(self):
        """ Bundles the writing/execution of simulation
        """
        self.write()
        self.partition()
        try:
            self.run_simulation()
        except subprocess.CalledProcessError as exc:
            print("Run failed, with output:")
            print(exc.output)
            raise
    
    def check(self):
        """ Load results and submit to check_results(dataset,grid)
        """
        map_ds=self.map_ds()
        g=unstructured_grid.UnstructuredGrid.read_ugrid(map_ds)
        self.check_results(map_ds,self.his_ds(),g)

    def check_results(self,map_ds,his_ds,g):
        print("No checks for %r"%self.__class__)

    def map_ds(self):
        map_fn=self.map_outputs()[0]
        ds=xr.open_dataset(map_fn)
        ds.close()
        return xr.open_dataset(map_fn)
    def his_ds(self):
        his_fn=self.his_output()
        his=xr.open_dataset(his_fn)
        his.close()
        return xr.open_dataset(his_fn)

    def plot(self):
        map_ds=self.map_ds()
        g=unstructured_grid.UnstructuredGrid.read_ugrid(map_ds)
        
        return self.plot_results(map_ds,self.his_ds(),g)

    def plot_results(self,map_ds,his_ds,g):
        fig=plt.figure(1)

        scalars=['mesh2d_ucmag']

        fig.clf()
        fig,axs=plt.subplots(len(scalars),1,sharex=True,num=self.fig_num,squeeze=False)
        axs=axs[:,0]

        tidx=-1
        for ax,scal in zip(axs,scalars):
            g.plot_edges(color='k',lw=0.5,ax=ax,)
            
            if scal=='age':
                values=map_ds['mesh2d_NO3'].isel(time=tidx).values / map_ds['mesh2d_RcNit'].isel(time=tidx).values
            else:
                values=map_ds[scal].isel(time=tidx)
            ccoll=g.plot_cells(values=values,cmap='jet',ax=ax)

            ax.axis('equal')
            plt.colorbar(ccoll,ax=ax,label=scal)
            
class AgeDomain(RectangleDomain):
    run_dir='run_age'
    fig_num=2
    
    def __init__(self,*a,**k):
        super(AgeDomain,self).__init__(*a,**k)

        # Basic age:
        self.dwaq.parameters["TcNit"]=1
        self.dwaq.parameters["NH4"]=1
        self.dwaq.substances['NO3']=self.dwaq.Sub(initial=0.0) # (name='NO3',active=True)
        self.dwaq.substances['RcNit']=self.dwaq.Sub(initial=0.0) # add_substance(name='RcNit',active=True)
        self.dwaq.add_process(name="Nitrif_NH4")
        
    def set_inflow(self):
        flow_bc=super(AgeDomain,self).set_inflow()
        conc_bc=dfm.DelwaqScalarBC(parent=flow_bc,scalar='RcNit',value=1)
        self.add_bcs([conc_bc])

    def check_results(self,map_ds,his_ds,g):
        ubar=map_ds['mesh2d_ucmag'].isel(time=-1).values.mean()
        tidx=-1
        age=map_ds['mesh2d_NO3'].isel(time=tidx).values / map_ds['mesh2d_RcNit'].isel(time=tidx).values
        # loose validation
        age_max_expected=(model.Lx) / ubar / 86400.
        rel_err=np.abs(age.max() - age_max_expected) / age_max_expected
        print("Age: max expected %.4f  max result %.4f"%( age_max_expected, age.max()))
        # Allow 5% error
        assert rel_err<0.05,"Age doesn't look right"
        
    def plot_results(self,map_ds,his_ds,g):
        scalars=['mesh2d_ucmag','mesh2d_RcNit','age']

        plt.figure(self.fig_num).clf()
        fig,axs=plt.subplots(len(scalars),1,sharex=True,num=self.fig_num,squeeze=False)
        axs=axs[:,0]

        tidx=-1
        for ax,scal in zip(axs,scalars):
            g.plot_edges(color='k',lw=0.5,ax=ax,)
            
            if scal=='age':
                values=map_ds['mesh2d_NO3'].isel(time=tidx).values / map_ds['mesh2d_RcNit'].isel(time=tidx).values
            else:
                values=map_ds[scal].isel(time=tidx)
            ccoll=g.plot_cells(values=values,cmap='jet',ax=ax)

            ax.axis('equal')
            plt.colorbar(ccoll,ax=ax,label=scal)

class PartialAgeDomain(AgeDomain):
    def check_results(self,map_ds,his_ds,g):
        # Skip AgeDomain
        super(AgeDomain,self).check_results(map_ds,his_ds,g)
        
        ubar=map_ds['mesh2d_ucmag'].isel(time=-1).values.mean()
        tidx=-1
        age=map_ds['mesh2d_NO3'].isel(time=tidx).values / map_ds['mesh2d_RcNit'].isel(time=tidx).values
        # loose validation
        age_max_expected=(model.Lx/2) / ubar / 86400.
        rel_err=np.abs(age.max() - age_max_expected) / age_max_expected
        print("Partial age: expected %.4f  results %.4f rel_err=%.3f"%(age_max_expected,age.max(),rel_err))
        # pretty generous bound here...
        assert rel_err<0.1,"Partial age doesn't look right"
        print("Max age looks okay")
        left_age=age[ g.cells_center()[:,0]<model.Lx/2 ]
        rel_err=np.abs(left_age/age_max_expected)
        assert rel_err.max()<0.01,"Shouldn't see age on the left"
        print("Left-size age looks okay")

class PartialPolyAgeDomain(PartialAgeDomain):
    """
    Tests specifying NH4 by polygon in order to get 
    partial age
    For partial age, use NH4 to select spatially
    This is supposed to end up in the ext file, such as waqparameterNH4
    Could have different types for value, and if it's a type that indicates
    spatial parameter, then it gets different handling in write_param() ?
    Manual says that a spatial or temporal variable set via .ext file
    overrules a constant added to the substances file.
    """
    fig_num=3
    run_dir="run_partial_poly"

    def write_forcing(self):
        super(PartialPolyAgeDomain,self).write_forcing()
        
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
        dio.write_pli(pli_fn, [ ('partial',np.array([ [500,0],
                                                      [1000,0],
                                                      [1000,50],
                                                      [500,50] ])) ])


class PartialRasterAgeDomain(PartialAgeDomain):
    """
    Tests specifying NH4 by constant in time raster
    """
    fig_num=4
    run_dir="run_partial_raster"
    
    def write_forcing(self):
        super(PartialRasterAgeDomain,self).write_forcing()
        
        bc_fn=self.ext_force_file()
        with open(bc_fn,'at') as fp:
            # filetype 11: netcdf gridded data
            # method 3: interpolate time and space, save coefficients
            # How to structure the nc file?
            # There is one example file in the distribution
            # "04_testcases/Example Testcases/f011_wind/c055_curvi_netcdf/wind_deltares_matroos_hirlam72_20150418.nc"

            stanza=f"""
QUANTITY=waqsegmentfunctionNH4
FILENAME=nh4.nc
VARNAME=nh4
FILETYPE=11
METHOD=3
OPERAND=O
"""
            fp.write(stanza)
        self.nh4_raster().to_netcdf( os.path.join(self.run_dir,'nh4.nc') )
            
    def nh4_raster(self):
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
        return ds

    
class SpatioTemporalAgeDomain(PartialRasterAgeDomain):
    """
    vary nh4 in time, too, by way of time-varying raster.
    """
    fig_num=5
    run_dir="run_spatiotemporal"

    def nh4_raster(self):
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
        return ds
    
    def check_results(self,map_ds,his_ds,g):
        super(SpatioTemporalAgeDomain,self).check_results(map_ds,his_ds,g)
        
        n_times=len(map_ds.time)
        t_idx=n_times//2 - 1
        age=map_ds['mesh2d_NO3'].isel(time=t_idx).values / map_ds['mesh2d_RcNit'].isel(time=t_idx).values

        assert age.max()<0.01,"Shouldn't see any age in first half of run"
        print("Time-varying age looks okay")


def test_rectangle():
    model=RectangleDomain(run_dir="run_rectangle")
    model.run_all()
    model.check()
    model.plot()

def test_age():
    model=AgeDomain(run_dir="run_age")
    model.run_all()
    model.check()
    model.plot()
    return model

def test_poly_age():
    model=PartialPolyAgeDomain()
    model.run_all()
    model.check()
    model.plot()
    return model

def test_raster_age():
    model=PartialRasterAgeDomain()
    model.run_all()
    model.check()
    model.plot()
    return model

def test_unsteady_raster_age():
    model=SpatioTemporalAgeDomain()
    model.run_all()
    model.check()
    model.plot()
    return model

##

# Run test_partial, but with a regenerated proc_def.def

class LocalAgeDomain(AgeDomain):
    run_dir='run_localage'
    fig_num=6

    proc_csv_dir="local_procs"

    def __init__(self,*a,**k):
        self.waq_proc_def=os.path.abspath( os.path.join(self.proc_csv_dir,"proc_def.def"))
        super(LocalAgeDomain,self).__init__(*a,**k)

        # Homegrown age:
        self.dwaq.substances['conc00']=self.dwaq.Sub()
        self.dwaq.substances['agec00']=self.dwaq.Sub()
        self.dwaq.add_process(name="AgeConc00")
    
    def set_inflow(self):
        flow_bc=super(AgeDomain,self).set_inflow() # skip AgeDomain implementation
        conc_bc=dfm.DelwaqScalarBC(parent=flow_bc,scalar='conc00',value=1)
        self.add_bcs([conc_bc])
    
    def run_simulation(self,threads=1,extra_args=[]):
        """
        Add in the Delwaq command line arguments
        """
        extra=extra_args+["--processlibrary",os.path.abspath(os.path.join(self.proc_csv_dir,"proc_def.def")),
                          "--openprocessdllso",os.path.abspath(os.path.join(self.proc_csv_dir,"liblocalprocess.so"))]
    
        return super(RectangleDomain,self).run_simulation(threads=threads,
                                                          extra_args=extra)

    def check_results(self,map_ds,his_ds,g):
        ubar=map_ds['mesh2d_ucmag'].isel(time=-1).values.mean()
        tidx=-1
        age=map_ds['mesh2d_agec00'].isel(time=tidx).values / map_ds['mesh2d_conc00'].isel(time=tidx).values
        # loose validation
        age_max_expected=(model.Lx) / ubar / 86400.
        rel_err=np.abs(age.max() - age_max_expected) / age_max_expected
        print("Age: max expected %.4f  max result %.4f"%( age_max_expected, age.max()))
        # Allow 5% error
        assert rel_err<0.05,"Age doesn't look right"
        
    def plot_results(self,map_ds,his_ds,g):
        scalars=['mesh2d_ucmag','mesh2d_conc00','age']

        plt.figure(self.fig_num).clf()
        fig,axs=plt.subplots(len(scalars),1,sharex=True,num=self.fig_num,squeeze=False)
        axs=axs[:,0]

        tidx=-1
        for ax,scal in zip(axs,scalars):
            g.plot_edges(color='k',lw=0.5,ax=ax,)
            
            if scal=='age':
                values=map_ds['mesh2d_agec00'].isel(time=tidx).values / map_ds['mesh2d_conc00'].isel(time=tidx).values
            else:
                values=map_ds[scal].isel(time=tidx)
            ccoll=g.plot_cells(values=values,cmap='jet',ax=ax)

            ax.axis('equal')
            plt.colorbar(ccoll,ax=ax,label=scal)

            
def test_local_age():
    model=LocalAgeDomain()
    model.run_all()
    model.check()
    model.plot()
    return model


#test_rectangle()
#test_age()
#test_poly_age()
#test_raster_age()
#test_unsteady_raster_age()
#test_local_age()
