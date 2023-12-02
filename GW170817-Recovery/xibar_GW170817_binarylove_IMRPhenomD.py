#!/usr/bin/env python
'''
Parameter estimation script that computes the posterior probability on all binary neutron star parameters, including xibar. 
Base waveform is IMRPhenomD.

The script makes use of the updated marginalized binary love relations arXiv:1903.03909.
'''
import sys
import bilby
import nrtidal_d
import numpy as np
from gwpy.timeseries import TimeSeries

#-----------------------------------------------------------------

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-o", "--outdir", type=str, default="outdir_GW170817_4",
                    help="Output directory for run.")

parser.add_argument("-l", "--label", type=str, default="bilby_GW170817_4",
                    help="Label of the bilby run.")

parser.add_argument("-n", "--npool", type=int, default=1, 
                    help="Number of CPUs.")

parser.add_argument("-sd", "--strain_dir", type=str, 
                    default = "/Users/abhi/Work/Projects/BDNK-Critical-Collapse/Data-Analysis/Cluster-Data/scratch/GW170817_files/Strain-Data-GW170817/no-glitch",
                    help="Path to directory that contains the strain data.")


args = parser.parse_args()

#-----------------------------------------------------------------

bilby.core.utils.setup_logger(outdir=args.outdir, label=args.label, log_level = "debug")
logger = bilby.core.utils.logger

#-----------------------------------------------------------------
'''
Set up sampling frequency and start, end times of the signal
'''
trigger_time = 1187008882.43

roll_off = 0.2  # Roll off duration of tukey window in seconds

# 4096 seconds of data

duration = 128   # Analysis segment duration
start_time = 1187008755
end_time = start_time + duration

sampling_frequency = 4096
#-----------------------------------------------------------------

'''
Read in the glitch free GW170817 data
'''

hdf5_filenames = {
    "H1": args.strain_dir+ "/H-H1_LOSC_CLN_4_V1-1187007040-2048_no_glitch.hdf5",
    "L1": args.strain_dir+ "/L-L1_LOSC_CLN_4_V1-1187007040-2048_no_glitch.hdf5",
    "V1": args.strain_dir+ "/V-V1_LOSC_CLN_4_V1-1187007040-2048_no_glitch.hdf5"
}

'''
Read in detector PSD files
'''

datapsd=np.loadtxt(args.strain_dir + '/GWTC1_GW170817_PSDs.dat')

farray=datapsd[:,0]
h1psd = datapsd[:,1]
l1psd= datapsd[:,2]
v1psd = datapsd[:,3]

psd_array = {"H1":h1psd, "L1":l1psd, "V1": v1psd}

det_names = np.array(["H1,L1,V1"])
ifo_list = bilby.gw.detector.InterferometerList([])

#Setting Strain Data
for det in ["H1", "L1", "V1"]:
    logger.info("Loading data for ifo {}".format(det))
    ifo = bilby.gw.detector.get_empty_interferometer(det)
    data = TimeSeries.read(
        hdf5_filenames[det], start=start_time,
        end=end_time, format="hdf5.gwosc"
    )
    ifo.strain_data.set_from_gwpy_timeseries(data)
    logger.info("Loading psd data for ifo {}".format(det))
    freq = farray
    psd = psd_array[det]
    logger.info("Setting IFO PSD from loaded PSD.")
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=freq, psd_array=psd)
    ifo_list.append(ifo)

logger.info("Finished setting up strain and PSD.")
logger.info("Saving IFO data plots to {}".format(args.outdir))
bilby.core.utils.check_directory_exists_and_if_not_mkdir(args.outdir)
ifo_list.plot_data(outdir=args.outdir, label=args.label)
#-------------------------------------------------------------
'''
Set up priors on all waveform parameters.
'''
priors = bilby.gw.prior.PriorDict()
priors['mass_1'] = bilby.gw.prior.Constraint(minimum=1, maximum=2,
                                             name='mass_1', latex_label='$m_1$', unit=None)
priors['mass_2'] = bilby.gw.prior.Constraint(minimum=1, maximum=2,
                                             name='mass_2', latex_label='$m_2$', unit=None)

priors['chirp_mass'] = bilby.gw.prior.UniformInComponentsChirpMass(minimum=1.184, maximum=1.25, name='chirp_mass', 
                                              latex_label='$\\mathcal{M}$', unit=None)
priors['mass_ratio'] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=0.5, maximum=1, name='mass_ratio',
                                              latex_label='$q$', unit=None)


priors['chi_1']=bilby.gw.prior.AlignedSpin(name='chi_1', a_prior=bilby.gw.prior.Uniform(minimum=0, maximum=0.05))
priors['chi_2']=bilby.gw.prior.AlignedSpin(name='chi_2', a_prior=bilby.gw.prior.Uniform(minimum=0, maximum=0.05))
priors['luminosity_distance'] = bilby.gw.prior.UniformSourceFrame(name='luminosity_distance',
                                                                     minimum=10, maximum=100, latex_label='$d_L$',
                                                                     unit='Mpc', boundary=None)
priors['phase'] = bilby.core.prior.Uniform(name='phase', minimum=0, maximum=2 * np.pi, boundary='periodic')
priors['theta_jn'] = bilby.prior.Sine(name='theta_jn', latex_label='$\\theta_{JN}$',
                                         unit=None, minimum=0, maximum=np.pi, boundary=None)
priors['psi'] = bilby.gw.prior.Uniform(name='psi', minimum=0, maximum=np.pi, boundary='periodic',
                                       latex_label='$\\psi$', unit=None)


priors["lambda_1"] = bilby.core.prior.Constraint( name="lambda_1", minimum=0, maximum=3000)
priors["lambda_2"] = bilby.core.prior.Constraint(name="lambda_2", minimum=0, maximum=3000)
priors["lambda_s"] = bilby.core.prior.Triangular(name="lambda_s", mode=1500,minimum = 0, maximum = 3000,latex_label='$\\Lambda_s$')

priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=trigger_time - 0.2,
    maximum=trigger_time + 0.2,
    name='geocent_time',
    latex_label='$t_c$',
    unit='$s$'
)
priors['dec'] =  bilby.prior.Cosine(name='dec', latex_label='$\\mathrm{DEC}$',
                                       unit=None, minimum=-np.pi / 2, maximum=np.pi / 2, boundary=None)
priors['ra'] =  bilby.gw.prior.Uniform(name='ra', minimum=0, maximum=2 * np.pi, boundary='periodic',
                                       latex_label='$\\mathrm{RA}$', unit=None)

priors["xi_bar"] = bilby.core.prior.Uniform(0,1000,name="xi_bar")

#-----------------------------------------------------------------
'''
Set up waveform generator and waveform source model
'''
waveform_arguments = dict(
    waveform_approximant="IMRPhenomD_NRTidal",
    reference_frequency=20.0
)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=nrtidal_d.source_binary_love,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    waveform_arguments=waveform_arguments
)

#-----------------------------------------------------------------
'''
Set up likelihood. We marginalize over the phase.
'''
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifo_list, waveform_generator=waveform_generator,
    time_marginalization=False, phase_marginalization=True,
    distance_marginalization=False, priors=priors)

#-----------------------------------------------------------------

result = bilby.run_sampler(
        likelihood=likelihood, 
        priors=priors, 
        sampler="dynesty", 
        sample = 'rwalk',
        bound = 'live',
        nlive=1500, 
        nact=10, 
        dlogz=0.01, 
        maxmcmc=5000,
        check_point_delta_t=3600,
        npool=args.npool, 
        outdir=args.outdir, label=args.label,
        conversion_function=bilby.gw.conversion.generate_all_bns_parameters)

result.plot_corner()

