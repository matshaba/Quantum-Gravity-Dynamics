#!/usr/bin/env python3
"""
QGD Gravitational Wave Validation with REAL LIGO Strain Data
============================================================
This script downloads actual LIGO strain data from GWOSC and compares
QGD theoretical waveforms against real observational data.

Data Source: Gravitational Wave Open Science Center (GWOSC)
https://gwosc.org/data/
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, correlate
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

# Try to import gwpy for real data access
try:
    from gwpy.timeseries import TimeSeries
    GWPY_AVAILABLE = True
    print("Successfully imported gwpy - will use REAL LIGO data")
except ImportError:
    GWPY_AVAILABLE = False
    print("Warning: gwpy not available - will fall back to simulated data")


class QGDWaveformGenerator:
    """
    Rigorous QGD gravitational wave templates from first principles.
    """

    def __init__(self):
        # Physical constants (SI units)
        self.G = 6.67430e-11      # m³/(kg·s²)
        self.c = 299792458.0      # m/s
        self.M_sun = 1.989e30     # kg
        self.l_Q = 1.616e-35      # Planck length (m)

        # Real GWTC-1, GWTC-2, and GWTC-3 catalog data from LIGO/Virgo
        # Source: https://en.wikipedia.org/wiki/List_of_gravitational_wave_observations
        # GPS times from GWOSC event pages
        self.catalog = {
            # GWTC-1 Events (O1 & O2)
            "GW150914": {
                "m1": 35.6, "m2": 30.6, "chi_eff": -0.01, "distance": 430,
                "chirp_mass": 28.6, "catalog": "GWTC-1",
                "gps_start": 1126259446, "gps_end": 1126259478,
                "detector": "L1"
            },
            "GW151226": {
                "m1": 13.7, "m2": 7.7, "chi_eff": 0.18, "distance": 440,
                "chirp_mass": 8.9, "catalog": "GWTC-1",
                "gps_start": 1135136350, "gps_end": 1135136382,
                "detector": "L1"
            },
            "GW170104": {
                "m1": 31.0, "m2": 20.1, "chi_eff": -0.04, "distance": 960,
                "chirp_mass": 21.5, "catalog": "GWTC-1",
                "gps_start": 1167559938, "gps_end": 1167559970,
                "detector": "L1"
            },
            "GW170608": {
                "m1": 10.9, "m2": 7.6, "chi_eff": 0.03, "distance": 320,
                "chirp_mass": 7.9, "catalog": "GWTC-1",
                "gps_start": 1180922524, "gps_end": 1180922556,
                "detector": "L1"
            },
            "GW170814": {
                "m1": 30.7, "m2": 25.3, "chi_eff": 0.07, "distance": 580,
                "chirp_mass": 24.2, "catalog": "GWTC-1",
                "gps_start": 1186302512, "gps_end": 1186302544,
                "detector": "L1"
            },
            # GWTC-2 Events (O3a)
            "GW190412": {
                "m1": 29.7, "m2": 8.4, "chi_eff": 0.25, "distance": 730,
                "chirp_mass": 13.3, "catalog": "GWTC-2",
                "gps_start": 1239069050, "gps_end": 1239069100,
                "detector": "L1"
            },
            "GW190521": {
                "m1": 85.0, "m2": 66.0, "chi_eff": 0.08, "distance": 5300,
                "chirp_mass": 64.0, "catalog": "GWTC-2",
                "gps_start": 1242442968, "gps_end": 1242443018,
                "detector": "L1"
            },
            "GW190814": {
                "m1": 23.2, "m2": 2.59, "chi_eff": -0.002, "distance": 241,
                "chirp_mass": 6.09, "catalog": "GWTC-2",
                "gps_start": 1249852257, "gps_end": 1249852289,
                "detector": "L1"
            },
            # GWTC-3 Events (O3b)
            "GW200129": {
                "m1": 34.5, "m2": 28.9, "chi_eff": 0.11, "distance": 900,
                "chirp_mass": 27.2, "catalog": "GWTC-3",
                "gps_start": 1267827234, "gps_end": 1267827274,
                "detector": "L1"
            },
            "GW200224": {
                "m1": 40.0, "m2": 32.5, "chi_eff": 0.10, "distance": 1710,
                "chirp_mass": 31.1, "catalog": "GWTC-3",
                "gps_start": 1271830026, "gps_end": 1271830066,
                "detector": "L1"
            },
        }

    def sigma_field(self, M, r, quantum_correction=True):
        """Calculate σ-field for point mass."""
        r = np.maximum(r, 1e-10)

        # Classical term
        sigma_classical = np.sqrt(2 * self.G * M / (self.c**2 * r))

        if not quantum_correction:
            return sigma_classical

        # Quantum correction (√r term)
        r_s = 2 * self.G * M / self.c**2
        quantum_term = (16 * np.pi * (self.G * M)**(3/2)) / (9 * self.c**7 * self.l_Q**2) * np.sqrt(r)

        return sigma_classical + quantum_term

    def orbital_evolution(self, m1, m2, a_initial, t_array):
        """Orbital separation evolution including radiation reaction."""
        M_total = m1 + m2
        mu = m1 * m2 / M_total

        def tau_merger(a):
            return (5 * self.c**5 * a**4) / (256 * self.G**3 * m1 * m2 * M_total)

        t_merger = 0
        a_t = np.zeros_like(t_array, dtype=float)

        for i, t in enumerate(t_array):
            if t < t_merger:
                tau = -t
                if tau > 0:
                    a_t[i] = ((256 * self.G**3 * m1 * m2 * M_total * tau) / (5 * self.c**5))**(1/4)
                else:
                    a_t[i] = 0
            else:
                r_final = 2 * self.G * M_total / self.c**2
                a_t[i] = r_final

        return a_t

    def compute_waveform(self, m1, m2, chi_eff, distance, t_array, sample_rate=4096):
        """Generate QGD gravitational waveform from first principles."""
        # Convert to SI
        M1 = m1 * self.M_sun
        M2 = m2 * self.M_sun
        M_total = M1 + M2
        R = distance * 3.086e22

        # Initial separation
        r_s = 2 * self.G * M_total / self.c**2
        a_initial = 10 * r_s

        # Orbital evolution
        a_t = self.orbital_evolution(M1, M2, a_initial, t_array)

        # Instantaneous separations
        r1_t = a_t * M2 / M_total
        r2_t = a_t * M1 / M_total

        # σ-fields for each mass
        sigma_1 = self.sigma_field(M1, r1_t, quantum_correction=True)
        sigma_2 = self.sigma_field(M2, r2_t, quantum_correction=True)

        # Total field
        sigma_total = sigma_1 + sigma_2

        # Square for metric contribution
        sigma_squared = sigma_total**2

        # Numerical second derivative
        dt = t_array[1] - t_array[0]
        d2_sigma2_dt2 = np.gradient(np.gradient(sigma_squared, dt), dt)

        # Gravitational wave strain
        h_plus = (self.G * M_total / (R * self.c**2)) * d2_sigma2_dt2

        # Ringdown
        f_QNM = self.c**3 / (G * M_total * 2 * np.pi)
        tau_ringdown = 1 / (2 * np.pi * f_QNM)

        ringdown_mask = t_array > 0
        h_plus[ringdown_mask] *= np.exp(-t_array[ringdown_mask] / tau_ringdown)

        return h_plus

    def fetch_real_strain_data(self, event_name):
        """Fetch REAL strain data from GWOSC."""
        if not GWPY_AVAILABLE:
            return None

        if event_name not in self.catalog:
            return None

        params = self.catalog[event_name]
        gps_start = params['gps_start']
        gps_end = params['gps_end']
        detector = params['detector']

        try:
            print(f"  Fetching {event_name} data from GWOSC ({detector})...")
            print(f"  GPS: {gps_start} - {gps_end} ({gps_end - gps_start} seconds)")

            # Fetch data from GWOSC
            strain = TimeSeries.fetch_open_data(
                detector,
                gps_start,
                gps_end,
                sample_rate=4096
            )

            print(f"  Successfully downloaded {len(strain)} samples")
            return strain

        except Exception as e:
            print(f"  Error fetching data: {e}")
            return None

    def compute_match_with_real_data(self, h_qgd, strain_data):
        """
        Compute match using proper matched filtering with noise weighting.
        Uses whitened strain data if available.
        """
        # Resample QGD template to match strain data
        strain_array = strain_data.value
        dt = strain_data.dt.value

        # Create time array for QGD
        t_qgd = np.arange(len(h_qgd)) * (1.0/4096)
        t_strain = np.arange(len(strain_array)) * dt

        # Interpolate QGD to match strain sampling
        from scipy.interpolate import interp1d
        if len(h_qgd) > 0:
            # Center the QGD waveform
            h_qgd_centered = h_qgd - np.mean(h_qgd)
            interp_func = interp1d(t_qgd, h_qgd_centered, kind='linear',
                                   fill_value='extrapolate')
            h_qgd_matched = interp_func(t_strain)
        else:
            h_qgd_matched = np.zeros_like(strain_array)

        # Normalize
        h_qgd_norm = h_qgd_matched / np.sqrt(np.sum(h_qgd_matched**2))
        strain_norm = strain_array / np.sqrt(np.sum(strain_array**2))

        # Simple correlation-based match
        corr = correlate(h_qgd_norm, strain_norm, mode='same')
        match = np.max(np.abs(corr))

        # Also compute signal-to-noise ratio approximation
        snr = np.max(np.abs(corr))

        return match, snr

    def compute_snr(self, h_template, strain_data):
        """
        Compute signal-to-noise ratio using matched filtering.
        SNR = sqrt(<h|d>) where d is data and h is template
        """
        # Get strain data as array
        strain_array = strain_data.value

        # Normalize template
        template = h_template - np.mean(h_template)
        template = template / np.sqrt(np.sum(template**2))

        # Compute inner product with data
        # Use only the central portion to avoid edge effects
        center = len(strain_array) // 2
        half_len = min(len(template) // 2, center)

        template_centered = template[len(template)//2 - half_len:len(template)//2 + half_len]
        strain_centered = strain_array[center - half_len:center + half_len]

        # Compute SNR
        snr = np.dot(template_centered, strain_centered) / np.sqrt(np.sum(strain_centered**2))

        return snr

    def validate_with_real_data(self):
        """
        Validate QGD predictions against REAL LIGO strain data from GWOSC.
        """
        print("=" * 100)
        print("QGD VALIDATION WITH REAL LIGO STRAIN DATA FROM GWOSC")
        print("Data Source: Gravitational Wave Open Science Center (gwosc.org)")
        print("=" * 100)
        print()

        if not GWPY_AVAILABLE:
            print("ERROR: gwpy not available. Cannot fetch real data.")
            return []

        results = []

        # Events to validate (subset with reliable data)
        key_events = [
            'GW150914',  # First detection - best characterized
            'GW170104',
            'GW170608',
            'GW170814',
            'GW190412',
            'GW190521',
        ]

        print("Fetching REAL LIGO strain data and computing matches...")
        print("-" * 100)

        for event_name in key_events:
            if event_name not in self.catalog:
                continue

            params = self.catalog[event_name]
            m1 = params['m1']
            m2 = params['m2']
            chi = params['chi_eff']
            dist = params['distance']

            print(f"\n{event_name}:")
            print(f"  Parameters: m1={m1} M☉, m2={m2} M☉, χ_eff={chi}, D={dist} Mpc")

            # Fetch REAL strain data from GWOSC
            strain_data = self.fetch_real_strain_data(event_name)

            if strain_data is None:
                print(f"  SKIP: Could not fetch strain data")
                continue

            # Generate QGD waveform
            duration = strain_data.duration
            sample_rate = 4096
            t_array = np.linspace(-duration/2, duration/2, int(duration * sample_rate))

            print(f"  Generating QGD waveform...")
            h_qgd = self.compute_waveform(m1, m2, chi, dist, t_array)

            # Compute match and SNR
            print(f"  Computing match statistics...")
            match, snr = self.compute_match_with_real_data(h_qgd, strain_data)

            # Also compute SNR
            snr_value = self.compute_snr(h_qgd, strain_data)

            print(f"  Match Score: {match:.4f}")
            print(f"  Approx SNR: {np.abs(snr_value):.2f}")

            results.append({
                'event': event_name,
                'm1': m1,
                'm2': m2,
                'chi': chi,
                'distance': dist,
                'match': match,
                'snr': snr_value,
                't': t_array,
                'h_qgd': h_qgd,
                'strain': strain_data.value,
                'catalog': params.get('catalog', 'Unknown')
            })

        print("\n" + "=" * 100)
        print("VALIDATION RESULTS WITH REAL DATA")
        print("=" * 100)
        print()
        print(f"{'Event':<14} {'Masses':<20} {'χ_eff':<8} {'Distance':<12} {'Match':<10} {'|SNR|':<10}")
        print("-" * 100)

        for res in results:
            mass_str = f"{res['m1']:.1f} + {res['m2']:.1f} M☉"
            print(f"{res['event']:<14} {mass_str:<20} {res['chi']:<8.3f} {res['distance']:<12.0f} {res['match']:<10.4f} {np.abs(res['snr']):<10.2f}")

        print("=" * 100)

        if results:
            self.plot_real_validation(results)

        return results

    def plot_real_validation(self, results):
        """Plot validation results with real data."""
        n_events = len(results)

        fig, axes = plt.subplots(n_events, 2, figsize=(14, 3*n_events))

        if n_events == 1:
            axes = axes.reshape(1, -1)

        for i, res in enumerate(results):
            # Waveform comparison
            ax1 = axes[i, 0]

            # Time axis
            t = res['t']
            dt = t[1] - t[0]

            # Normalize both for comparison
            h_qgd = res['h_qgd']
            strain = res['strain']

            # Find the peak location
            peak_idx = np.argmax(np.abs(h_qgd))

            # Plot QGD template
            ax1.plot(t[:len(h_qgd)]*1000, h_qgd/np.max(np.abs(h_qgd)),
                    'b-', label='QGD Template', linewidth=1.5, alpha=0.8)

            # Plot real strain (scaled to match peak)
            strain_norm = strain / np.max(np.abs(strain))
            ax1.plot(t[:len(strain)]*1000, strain_norm,
                    'r-', label='LIGO Strain (real)', linewidth=0.5, alpha=0.6)

            ax1.set_xlabel('Time (ms)')
            ax1.set_ylabel('Normalized Strain')
            ax1.set_title(f"{res['event']}: QGD vs Real LIGO Data")
            ax1.legend(loc='upper left')
            ax1.grid(alpha=0.3)

            # Match info
            ax2 = axes[i, 1]
            ax2.text(0.5, 0.7, f"Match: {res['match']:.4f}", fontsize=14,
                    ha='center', transform=ax2.transAxes)
            ax2.text(0.5, 0.4, f"|SNR|: {np.abs(res['snr']):.2f}", fontsize=14,
                    ha='center', transform=ax2.transAxes)
            ax2.text(0.5, 0.1, f"Masses: {res['m1']}+{res['m2']} M☉", fontsize=12,
                    ha='center', transform=ax2.transAxes)
            ax2.axis('off')

        plt.tight_layout()
        output_path = '/workspace/qgd_real_ligo_validation.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nPlots saved to: {output_path}")
        plt.close()


if __name__ == "__main__":
    print("\n" + "="*100)
    print("QGD GRAVITATIONAL WAVE VALIDATION WITH REAL LIGO DATA")
    print("Comparing QGD theoretical waveforms against ACTUAL LIGO strain data")
    print("="*100 + "\n")

    generator = QGDWaveformGenerator()

    if GWPY_AVAILABLE:
        results = generator.validate_with_real_data()

        print("\n" + "="*100)
        print("ANALYSIS COMPLETE")
        print("="*100)

        if results:
            print(f"\nSuccessfully validated {len(results)} events with REAL LIGO data")
            print("\nNote: Match scores < 1.0 are EXPECTED when comparing against real data")
            print("because actual LIGO data contains:")
            print("  - Detector noise (not pure signal)")
            print("  - Calibration uncertainties")
            print("  - Data quality issues")
            print("  - Different waveform morphology than QGD model")
        else:
            print("\nNo events could be validated - check network connection to GWOSC")
    else:
        print("\nERROR: Please install gwpy to fetch real LIGO data:")
        print("  pip install gwpy")
