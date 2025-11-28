"""
Gravitational Wave Dispersion Simulation
for the E↔S Framework Paper

Simulates the propagation of a GW pulse with massive graviton dispersion relation:
    E² = p²c² + m_g²c⁴
    
This leads to frequency-dependent group velocity:
    v_g(f) = c * sqrt(1 - (m_g c² / h f)²)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, h, eV

# Physical constants
c_light = c  # m/s
h_planck = h  # J·s
eV_to_J = eV  # J per eV

# Source parameters (typical binary merger)
D_source = 400e6 * 3.086e16  # 400 Mpc in meters
f_min = 20  # Hz (LIGO lower bound)
f_max = 500  # Hz (typical merger frequency)

def group_velocity(f, m_g_eV):
    """
    Calculate group velocity for gravitational waves with massive graviton.
    
    Parameters:
    -----------
    f : float or array
        Frequency in Hz
    m_g_eV : float
        Graviton mass in eV/c² (so m_g c² is in eV)
    
    Returns:
    --------
    v_g : group velocity in m/s
    """
    E = h_planck * f  # Energy in Joules
    m_g_c2_J = m_g_eV * eV_to_J  # m_g c² in Joules (rest mass energy)
    
    # v_g = c * sqrt(1 - (m_g c² / E)²)
    ratio_squared = (m_g_c2_J / E)**2
    
    # Avoid sqrt of negative (frequencies below threshold)
    ratio_squared = np.clip(ratio_squared, 0, 0.9999)
    
    return c_light * np.sqrt(1 - ratio_squared)

def time_delay(f, m_g_eV, D):
    """
    Calculate arrival time delay relative to massless case.
    
    For m_g c² << E: Δt ≈ D/(2c) * (m_g c² / E)²
    
    This approximation is numerically stable for small mass ratios.
    """
    E = h_planck * f  # Energy in Joules
    m_g_c2_J = m_g_eV * eV_to_J  # m_g c² in Joules
    
    # Use approximation: Δt ≈ D/(2c) * (m_g c² / E)²
    ratio = m_g_c2_J / E
    return (D / (2 * c_light)) * ratio**2

def simulate_waveform(t, f_chirp, phase=0):
    """Simple chirp waveform for visualization."""
    return np.sin(2 * np.pi * f_chirp * t + phase) * np.exp(-((t - t.mean())/0.05)**2)

# =============================================================================
# SIMULATION 1: Time delay vs frequency for different graviton masses
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Mass values to explore (LIGO limit is ~10^-22 eV)
mass_values = [1e-23, 1e-22, 1e-21, 1e-20]  # eV
frequencies = np.linspace(f_min, f_max, 1000)

ax1 = axes[0, 0]
for m_g in mass_values:
    delays = time_delay(frequencies, m_g, D_source)
    ax1.plot(frequencies, delays, label=f'$m_g = 10^{{{int(np.log10(m_g))}}}$ eV')

ax1.set_xlabel('Frequency (Hz)', fontsize=12)
ax1.set_ylabel('Time delay Δt (seconds)', fontsize=12)
ax1.set_title('GW Arrival Time Delay vs Frequency\n(Source at 400 Mpc)', fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# =============================================================================
# SIMULATION 2: Group velocity ratio v_g/c (using approximation for small masses)
# =============================================================================

ax2 = axes[0, 1]
for m_g in mass_values:
    # For small masses: 1 - v_g/c ≈ (1/2)(m_g c² / E)²
    E = h_planck * frequencies
    m_g_c2_J = m_g * eV_to_J
    v_reduction = 0.5 * (m_g_c2_J / E)**2
    ax2.plot(frequencies, v_reduction, label=f'$m_g = 10^{{{int(np.log10(m_g))}}}$ eV')

ax2.set_xlabel('Frequency (Hz)', fontsize=12)
ax2.set_ylabel('$(c - v_g)/c$', fontsize=12)
ax2.set_title('Fractional Velocity Reduction', fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

# =============================================================================
# SIMULATION 3: Waveform dispersion visualization
# =============================================================================

ax3 = axes[1, 0]

# Create a GW pulse with two frequency components
t = np.linspace(-0.2, 0.2, 2000)
f_low = 50  # Hz
f_high = 200  # Hz

m_g_demo = 1e-22  # eV (at LIGO limit for visible effect)

# Time delays for each component
dt_low = time_delay(f_low, m_g_demo, D_source)
dt_high = time_delay(f_high, m_g_demo, D_source)
delta_t = dt_low - dt_high  # Relative delay

# Massless case (both arrive together)
signal_massless = (0.7 * simulate_waveform(t, f_low) + 
                   0.5 * simulate_waveform(t, f_high))

# Massive case (low frequency delayed)
signal_massive = (0.7 * simulate_waveform(t - delta_t/2, f_low) + 
                  0.5 * simulate_waveform(t + delta_t/2, f_high))

ax3.plot(t * 1000, signal_massless, 'b-', alpha=0.7, label='Massless graviton')
ax3.plot(t * 1000, signal_massive, 'r--', alpha=0.7, label=f'$m_g = 10^{{-22}}$ eV (LIGO limit)')
ax3.set_xlabel('Time (ms)', fontsize=12)
ax3.set_ylabel('Strain h(t)', fontsize=12)
ax3.set_title(f'Waveform Dispersion Effect\nΔt = {delta_t*1000:.2f} ms between {f_low} Hz and {f_high} Hz', fontsize=12)
ax3.legend()
ax3.grid(True, alpha=0.3)

# =============================================================================
# SIMULATION 4: Constraint plot - Observable delay vs mass
# =============================================================================

ax4 = axes[1, 1]

masses = np.logspace(-24, -19, 100)  # eV
delays_50Hz = [time_delay(50, m, D_source) for m in masses]
delays_100Hz = [time_delay(100, m, D_source) for m in masses]

ax4.loglog(masses, delays_50Hz, 'b-', label='f = 50 Hz')
ax4.loglog(masses, delays_100Hz, 'r-', label='f = 100 Hz')
ax4.axvline(x=1e-22, color='k', linestyle='--', alpha=0.5, label='LIGO limit')
ax4.axhline(y=1e-3, color='gray', linestyle=':', alpha=0.5, label='1 ms threshold')

ax4.set_xlabel('Graviton mass $m_g$ (eV)', fontsize=12)
ax4.set_ylabel('Time delay (seconds)', fontsize=12)
ax4.set_title('Observable Time Delay vs Graviton Mass\n(Source at 400 Mpc)', fontsize=12)
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_xlim(1e-24, 1e-19)
ax4.set_ylim(1e-6, 1e6)

plt.tight_layout()
plt.savefig('gw_dispersion_simulation.png', dpi=150, bbox_inches='tight')
plt.savefig('gw_dispersion_simulation.pdf', bbox_inches='tight')
print("Figures saved: gw_dispersion_simulation.png/pdf")

# =============================================================================
# NUMERICAL RESULTS
# =============================================================================

print("\n" + "="*60)
print("NUMERICAL RESULTS FOR PAPER")
print("="*60)

print(f"\nSource distance: D = 400 Mpc = {D_source:.2e} m")
print(f"LIGO frequency range: {f_min} - {f_max} Hz")

print("\n--- Time delays at f = 100 Hz ---")
for m_g in [1e-23, 1e-22, 1e-21]:
    dt = time_delay(100, m_g, D_source)
    print(f"m_g = 10^{int(np.log10(m_g)):3d} eV:  Δt = {dt:.3e} s = {dt*1000:.3f} ms")

print("\n--- Differential delay between 50 Hz and 200 Hz ---")
for m_g in [1e-23, 1e-22, 1e-21]:
    dt_diff = time_delay(50, m_g, D_source) - time_delay(200, m_g, D_source)
    print(f"m_g = 10^{int(np.log10(m_g)):3d} eV:  Δt(50Hz) - Δt(200Hz) = {dt_diff:.3e} s")

print("\n--- Critical note for E↔S paper ---")
print("The paper claims m_Φ = M_P ~ 10^19 GeV = 10^28 eV")
print("This is 50 orders of magnitude above LIGO limits!")
print("For m_g ~ M_P, GWs with f < 10^43 Hz would not propagate.")
print("The screening mechanism must reduce the effective mass dramatically.")

plt.show()
