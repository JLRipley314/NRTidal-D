#!/usr/bin/env python

import sys
import bilby
import nrtidal_d
import numpy as np

#-----------------------------------------------------------------

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-o", "--outdir", type=str, default="output",
                    help="Output directory for run.")

parser.add_argument("-l", "--label", type=str, default="GW170817",
                    help="Label of the bilby run.")

parser.add_argument("-n", "--npool", type=int, default=1, 
                    help="Number of CPUs.")

parser.add_argument("-x", "--xitilde", type=float, default=0, 
                    help="Value of xi tilde")

parser.add_argument("-ad", "--asd_dir", type=str, default = "", 
                    help="Location of strain data")

args = parser.parse_args()

#-----------------------------------------------------------------

bilby.core.utils.setup_logger(outdir=args.outdir, label=args.label)#, log_level = "debug")
logger = bilby.core.utils.logger

#-----------------------------------------------------------------
trigger_time = 1187008882.43

roll_off = 0.2  # Roll off duration of tukey window in seconds

duration = 128   # Analysis segment duration
start_time = 1187008755
end_time = start_time + duration

sampling_frequency = 2048 
#-----------------------------------------------------------------
#
# GW170817-like event
#
injection_parameters = dict(
        mass_1=1.38,
        mass_2=1.38,
        a_1=0.0,
        a_2=0.0,
        tilt_1=0.0,
        tilt_2=0.0,
        luminosity_distance=40.0,
        geocent_time=1187008882.4,
        theta_jn=2.64,
        psi=1.8,
        phase=0.0,
        ra=3.4,
        dec=-0.401,
        lambda_s=584,
        xi_tilde=args.xitilde,
        fiducial=1,
    )
#-----------------------------------------------------------------
priors = bilby.core.prior.PriorDict()
for key in list(injection_parameters.keys()):
    priors[key] = injection_parameters[key]
    
del priors["mass_1"], priors["mass_2"]

priors["chirp_mass"] = bilby.gw.prior.UniformInComponentsChirpMass(
    1.18, 2.17, name="chirp_mass", unit="$M_{\\odot}$"
)
priors["mass_ratio"] = bilby.gw.prior.UniformInComponentsMassRatio(
    0.1, 1.0, name="mass_ratio"
)

priors["lambda_1"] = bilby.core.prior.Constraint(
    name="lambda_1", minimum=0, maximum=3000
)
priors["lambda_2"] = bilby.core.prior.Constraint(
    name="lambda_2", minimum=0, maximum=3000
)
priors["lambda_s"] = bilby.core.prior.Triangular(mode=1500,minimum=0,maximum=3000)

priors["xi_tilde"] = bilby.core.prior.Uniform(0,1000,name="xi_tilde")

del priors["geocent_time"], priors["phase"]
#-----------------------------------------------------------------
waveform_arguments = dict(
    waveform_approximant="IMRPhenomPv2_NRTidal",
    reference_frequency=50.0,
    minimum_frequency=40.0
)
#-----------------------------------------------------------------
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=nrtidal_d.source_binary_love,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    waveform_arguments=waveform_arguments
)
#-----------------------------------------------------------------
CE_asd = np.loadtxt(args.asd_dir + "cosmic_explorer_strain.txt")

CE_f = CE_asd[:,0]

CE_asd = CE_asd[:,1]

f_array   = {"CE":CE_f}
asd_array = {"CE":CE_asd}

ifo_list = bilby.gw.detector.InterferometerList([])

#Setting Strain Data
for det in ["CE"]:
    logger.info("Loading asd data for ifo {}".format(det))
    ifo = bilby.gw.detector.get_empty_interferometer(det)
    
    freq = f_array[det]  
    asd = asd_array[det]
    
    logger.info("Setting IFO ASD from loaded ASD.")
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=freq, asd_array=asd)
    ifo_list.append(ifo)

for ifo in ifo_list:
    ifo.minimum_frequency = waveform_arguments["minimum_frequency"]

ifo_list.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency, duration=duration, start_time=start_time
)
ifo_list.inject_signal(
    parameters=injection_parameters, waveform_generator=waveform_generator
)

logger.info("Finished setting up strain and ASD.")
logger.info("Saving IFO data plots to {}".format(args.outdir))
bilby.core.utils.check_directory_exists_and_if_not_mkdir(args.outdir)
ifo_list.plot_data(outdir=args.outdir, label=args.label)

#-----------------------------------------------------------------
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifo_list, 
    waveform_generator=waveform_generator,
    time_marginalization=True, 
    phase_marginalization=True,
    distance_marginalization=False, 
    priors=priors)

#-----------------------------------------------------------------
result = bilby.run_sampler(
        likelihood=likelihood, 
        priors=priors, 
        sampler="dynesty", 
        sample = 'rwalk',
        bound = 'live',
        nlive=1500, 
        nact=5, 
        dlogz=0.1, 
        maxmcmc=5000,
        check_point_delta_t=7200,
        npool=args.npool, 
        outdir=args.outdir, label=args.label,
        conversion_function=bilby.gw.conversion.generate_all_bns_parameters)

result.plot_corner()
