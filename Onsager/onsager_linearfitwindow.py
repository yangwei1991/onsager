import os
import numpy as np
import pandas as pd
from scipy.stats import linregress
import mdtraj as mt
import matplotlib.pyplot as plt

def convert_to_gro(file, outfile, cell):
    '''
    Convert .dat file to .gro so it can be read by mdtraj; mdtraj calculates delta r vectors accounting for pbc
    
    file (str): .dat file
    outfile (str): .gro file
    cell (list of float): unit cell lengths in nm
    '''
    data = np.loadtxt(file, comments=['#'])
    xyz_list = data[:, :]
    with open(outfile, 'w+') as f:
        f.write('Comment\n')
        f.write(f'{len(xyz_list)}\n')
        for i, xyz in enumerate(xyz_list):
            x, y, z = xyz/10  # Convert Å to nm
            f.write(f'{i+1:>5}H        H{i+1:>5}{x:>8.3f}{y:>8.3f}{z:>8.3f}\n')
        f.write(f'{cell[0]} {cell[1]} {cell[2]}\n')  # Set pbc in .gro
        
def convert_to_gro_all(ddir, prefix, n_molecules, conc, cell):
    '''
    Convert all .dat files to .gro
    
    ddir (str): directory containing .dat files for all ions
    prefix (str): ion name for .dat files
    n_molecules (int): number of ions; each ion has a .dat file
    conc (str): concentration indicator for .dat files
    cell (list of float): unit cell lengths in nm
    '''
    for i in range(n_molecules):
        file = os.path.join(ddir, f'{prefix}{i+1}_{conc}.dat')
        outfile = f'{prefix}{i+1}.gro'
        convert_to_gro(file, outfile, cell)

def get_dr(traj, dt):
    '''
    Calculate delta r (vector) at all t for a given delta t and a given molecule
    
    traj (mt.Trajectory): mdtraj object
    dt (int): delta t in saved timesteps
    '''
    last_t = traj.n_atoms - dt
    # Get all timesteps
    indices = [[t, t+dt] for t in range(last_t)]
    # Calculate delta r for given delta at all timesteps
    dr = mt.compute_displacements(traj, indices, periodic=True)
    return dr.reshape(-1, 3)

def get_dr_all(prefix, n_molecules, dt_max):
    '''
    Calculate delta r at all t for a given range of delta t; then sums across all molecules of same type
    
    prefix (str): name indicating ion
    n_molecules (int): number of ions; each ion has a .dat file
    dt_max (int): maximum delta t to calculate (in saved timesteps)
    conc (str): concentration indicator for .dat files
    '''
    # Calculate delta r at all t for delta t between 1 and dt_max (saved timesteps)
    # Delta r calculated for all molecules of the same type
    dr_all = [] # L x M x N x 3 array of dr values; L is number of molecules; M is number of delta t; N is number of t
    outfiles = []
    for i in range(n_molecules):
        outfile = f'{prefix}{i+1}.gro'
        traj = mt.load(outfile)
        dr_list = [get_dr(traj, dt) for dt in range(1, dt_max)]
        dr_all.append(dr_list)
        outfiles.append(outfile)
    # Sum delta r at each delta t for all molecules of the same type
    dr_sum = [] # M x N x 3 array of dr values; M is number of delta t; N is number of t
    for j in range(dt_max-1):
        dr_sum.append(np.sum([dr_all[i][j] for i in range(n_molecules)], axis=0))
    # Clean up .gro files
    for outfile in outfiles:
        os.system(f'rm {outfile}')
    return dr_sum

def get_br(dr_sum1, dr_sum2):
    '''
    Calculate bracket term for each delta t by averaging over t
    
    dr_sum1 (list of np.array): M x N x 3 array of dr values; M is number of delta t; N is number of t
    dr_sum2 (list of np.array): M x N x 3 array of dr values; M is number of delta t; N is number of t
    '''
    # Calculate inner product and then average across all t for each delta t
    br = []
    for dr1, dr2 in zip(dr_sum1, dr_sum2):
        prod = np.sum(dr1*dr2, axis=-1)
        br.append(np.mean(prod))
    return br

def get_L(time, br, dt_end, symbol, cell, T):
    '''
    Get Onsager coefficient (L), e^2*L (s), and correlation coefficient (r)
    
    time (list): list of delta_t values (in ps)
    br (list): list of bracket terms
    dt_end (int): maximum delta t to perform linear fit (in saved timesteps)
    symbol (str): ++, --, or +- to show in linear fit plot
    cell (list of float): unit cell lengths in nm
    T (float): temperature in K
    '''
    # Linear fit between bracket term and delta t
    result = linregress(time[:dt_end], br[:dt_end])
    slope = result[0]  # units: nm^2 / ps
    intercept = result[1]
    r = result[2]
    # Get fitted line
    br_hat = slope*time[:dt_end] + intercept
    # Calculate L
    kB = 1.380649e-23
    V = cell[0]*cell[1]*cell[2] # units: nm^3
    L = 1/(6*kB*T*V)*slope # units: J-1 nm-1 ps-1
    # Convert to SI units
    L = L*1/(1e-9)*1/(1e-12) # units: J-1 m-1 s-1
    # Calculate e^2*L
    e = 1.602176634e-19
    s = e**2 * L # units: S/m
    # Plot bracket term vs. delta t along with linear fit
    plot_L(time[:dt_end], br[:dt_end], br_hat, slope, r, L, s, symbol)
    return r, L, s

def plot_L(time, br, br_hat, slope, r, L, s, symbol):
    '''
    Plot bracket term vs. delta t along with linear fit
    
    time (list): list of delta_t values (in ps)
    br (list): list of bracket terms
    br_hat (list): list of bracket terms calculated from linear fit
    slope (float): slope of linear fit
    L (float): Onsager coefficient
    s (float): e^2*L
    r (float): correlation coefficient
    symbol (str): ++, --, or +- to show in linear fit plot
    '''
    fig, ax = plt.subplots(1, 1)
    ax.plot(time, br, color='navy', linewidth=2, alpha=0.5)
    ax.plot(time, br_hat, '--', color='black', linewidth=2)
    # Show slope value on plot
    ax.text(0.5, 0.4, f'slope = {slope:.4f}'r' nm$^2$ ps$^{-1}$', transform=ax.transAxes)
    # Show r value on plot
    ax.text(0.5, 0.3, r'$r$ = 'f'{r:.3f}', transform=ax.transAxes)
    # Show L value on plot
    ax.text(0.5, 0.2, r'$L^{REP}$ = '.replace('REP', symbol) + f'{L:.3e}'r' J$^{-1}$ m$^{-1}$ s$^{-1}$', 
            transform=ax.transAxes)
    # Show e^2*L value on plot
    ax.text(0.5, 0.1, r'$e^2L^{REP}$ = 'f'{s:.3e}'.replace('REP', symbol) + r' S/m', transform=ax.transAxes)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel(r'Bracket term (nm$^2$)')
    ax.set_title(r' $L^{REP}$ fitting'.replace('REP', symbol))
    fig.tight_layout()
    fig.set_dpi(100)
    
def get_onsager_modified(ddir, n1=4, n2=8, dt_max_init=3, dt_time=50., prefix1='CA', prefix2='TFS', 
                         conc='0.05M', z1=2., z2=-1., cell=[5.2507, 5.2507, 5.2507], T=300., r_min=0.99, r_max=0.999):
    '''
    Adjust dt_max and dt_end such that 0.99 < r < 0.999 for correlation coefficients
    
    ddir (str): directory with .dat files
    n1 (int): number of cations
    n2 (int): number of anions
    dt_max_init (int): initial maximum delta t to calculate (in saved timesteps)
    dt_time (float): time in ps to which the saved timestep (delta t interval) corresponds (default is 50 ps)
    prefix1 (str): cation name for .dat files
    prefix2 (str): anion name for .dat files
    conc (str): concentration indicator for .dat files
    z1 (float): charge of cation
    z2 (float): charge of anion (including the negative sign, e.g., z2 = -1)
    cell (list of float): unit cell lengths in nm
    T (float): temperature in K
    r_min (float): minimum acceptable correlation coefficient
    r_max (float): maximum acceptable correlation coefficient
    '''
    dt_max = dt_max_init
    while dt_max <= 100:  # Loop to adjust dt_max and dt_end values up to 100 timesteps
        convert_to_gro_all(ddir, prefix=prefix1, n_molecules=n1, conc=conc, cell=cell)
        convert_to_gro_all(ddir, prefix=prefix2, n_molecules=n2, conc=conc, cell=cell)

        dr_sum1 = get_dr_all(prefix=prefix1, n_molecules=n1, dt_max=dt_max)
        dr_sum2 = get_dr_all(prefix=prefix2, n_molecules=n2, dt_max=dt_max)
        
        br_ii = get_br(dr_sum1, dr_sum1)
        br_jj = get_br(dr_sum2, dr_sum2)
        br_ij = get_br(dr_sum1, dr_sum2)

        time = np.arange(1, dt_max) * dt_time

        rii, Lii, sii = get_L(time, br_ii, dt_end=dt_max, symbol='++', cell=cell, T=T)
        rjj, Ljj, sjj = get_L(time, br_jj, dt_end=dt_max, symbol='--', cell=cell, T=T)
        rij, Lij, sij = get_L(time, br_ij, dt_end=dt_max, symbol='+-', cell=cell, T=T)

        if all(r_min < r < r_max for r in [rii, rjj, rij]):
            sigma = (z1**2 * sii + z2**2 * sjj + 2 * z1 * z2 * sij) * 10  # Conductivity in mS/cm
            return sigma, ((rii, Lii, sii), (rjj, Ljj, sjj), (rij, Lij, sij)), dt_max

        dt_max += 1

    print(f"No suitable dt_max found within range. Last values used: dt_max={dt_max_init}")
    return None

if __name__ == "__main__":
    sigma, data, dt_max_used = get_onsager_modified(ddir='data_files', n1=4, n2=8, dt_max_init=3, dt_time=50., 
                                                    prefix1='CA', prefix2='TFS', conc='0.05M', z1=2., z2=-1., 
                                                    cell=[5.2507]*3, T=300., r_min=0.99, r_max=0.999)

if sigma:
    with open(f"output_results_{dt_max_used}.txt", "w") as f:
            f.write(f"Conductivity: {sigma} mS/cm\n")
            f.write(f"Optimal dt_max: {dt_max_used}\n")
            
            r_ii, L_ii, s_ii = data[0]
            r_jj, L_jj, s_jj = data[1]
            r_ij, L_ij, s_ij = data[2]
            
            f.write(f"\nCations (CA):\n")
            f.write(f"Correlation coefficient: {r_ii}\n")
            f.write(f"Value of e^2 * L: {s_ii} S/m\n")
            
            f.write(f"\nAnions (TFS):\n")
            f.write(f"Correlation coefficient: {r_jj}\n")
            f.write(f"Value of e^2 * L: {s_jj} S/m\n")
            
            f.write(f"\nCross-term (Cations-Anions):\n")
            f.write(f"Correlation coefficient: {r_ij}\n")
            f.write(f"Value of e^2 * L: {s_ij} S/m\n")