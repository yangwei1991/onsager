import os
import numpy as np
from scipy.optimize import minimize
import json
import subprocess

# Function to run energy minimization and return energy
def run_minimization(params, folder, structure_count=27):
    l, n, p, r, t, v, c1 = params
    energies = []

    os.chdir(folder)  # Change directory
    print(f"Changed directory to: {os.getcwd()}")  # Logging current directory

    prmtop_filename = f"s000001_dry_{l:.4f}_{n:.4f}_{p:.4f}_{r:.4f}_{t:.4f}_{v:.4f}_{c1:.4f}.prmtop"

    # Modify LJ parameters using parmed on s000001_dry.prmtop
    parmed_input = f"""
    changeLJSingleType @1 1.3600 {l}
    changeLJSingleType @2 1.8993 {n}
    changeLJSingleType @3,10 1.9069 {p}
    changeLJSingleType @4,7 2.2777 {r}
    changeLJSingleType @5,6,8,9 1.7107 {t}
    changeLJSingleType @11-16 1.7029 {v}
    change CHARGE @1 {2.0000+c1}
    change CHARGE @2 {-0.466327+c1*-0.466327}
    change CHARGE @3 {0.331932+c1*0.331932}
    change CHARGE @4 {0.569588+c1*0.569588}
    change CHARGE @5 {-0.403534+c1*-0.403534}
    change CHARGE @6 {-0.372787+c1*-0.372787}
    change CHARGE @7 {0.569698+c1*0.569698}
    change CHARGE @8 {-0.372762+c1*-0.372762}
    change CHARGE @9 {-0.403258+c1*-0.403258}
    change CHARGE @10 {0.331854+c1*0.331854}
    change CHARGE @11 {-0.113780+c1*-0.113780}
    change CHARGE @12 {-0.166440+c1*-0.166440}
    change CHARGE @13 {-0.112236+c1*-0.112236}
    change CHARGE @14 {-0.165822+c1*-0.165822}
    change CHARGE @15 {-0.113850+c1*-0.113850}
    change CHARGE @16 {-0.112279+c1*-0.112279}
    outparm {prmtop_filename}
    quit
    """

    with open('parmed.in', 'w') as f:
        f.write(parmed_input)

    # Run parmed to modify parameters
    os.system(f"parmed -p s000001_dry.prmtop -i parmed.in")

    for a in range(1, structure_count + 1):
        min_out_filename = f"min_{a:02d}_{l:.4f}_{n:.4f}_{p:.4f}_{r:.4f}_{t:.4f}_{v:.4f}_{c1:.4f}.out"

        # Run energy minimization using sander
        os.system(f"sander -O -i min.in -o {min_out_filename} -p {prmtop_filename} -c s0000{a:02d}_dry.inpcrd -r min_{a:02d}.rst -x min_{a:02d}.netcdf")

        # Extract energy using grep, tail, head, and awk
        energy = None
        try:
            grep_command = f"grep -A 5 'FINAL RESULTS' {min_out_filename} | tail -n +6 | head -n 1 | awk '{{print $2}}'"
            energy_output = subprocess.check_output(grep_command, shell=True).decode().strip()
            energy = float(energy_output)
            energies.append(energy)
        except (subprocess.CalledProcessError, ValueError) as e:
            print(f"{e}")

        # Cleanup unnecessary files
        cleanup_files = [
            f"min_{a:02d}.rst",
            f"min_{a:02d}.netcdf",
                                 min_out_filename,
        ]
        for file in cleanup_files:
            if os.path.isfile(file):
                os.remove(file)

        if energy is None:
            raise RuntimeError(f"Failed to extract energy for params: {l, n, p, r, t, v, c1}")

        print(f"Energy for {os.path.basename(folder)} structure {a}: {energy}")

    # Remove the modified prmtop file
    if os.path.isfile(prmtop_filename):
        os.remove(prmtop_filename)
        print(f"Removed modified prmtop file: {prmtop_filename}")

    os.chdir('..')  # Move back to the parent directory
    print(f"Returned to parent directory: {os.getcwd()}")  # Logging current directory

    return energies

# Function to read QM energies from the file
def read_qm_energies():
    with open('QM_each_ComplexEnergy_relativeenergy.txt', 'r') as f:
        qm_energies = [float(line.strip()) for line in f.readlines()]
    return qm_energies

# Objective function for optimization
def objective_function(params, folders):
    qm_energies = read_qm_energies()
    mm_energies = []

    for folder in folders:
        energies = run_minimization(params, folder)
        mm_energies.extend(energies)

    # Use the lowest energy as the reference
    ref_energy = min(mm_energies)
    relative_mm_energies = [energy - ref_energy for energy in mm_energies]

#    # Save all energies to a file
#    energy_filename = f"energy_{params[0]:.4f}_{params[1]:.4f}_{params[2]:.4f}_{params[3]:.4f}_{params[4]:.4f}_{params[5]:.4f}_{params[6]:.4f}_{params[7]:.4f}.txt"
#    with open(energy_filename, 'w') as f:
#        for i, energy in enumerate(relative_mm_energies):
#            f.write(f"{energy}\n")

    # Calculate the difference between MM and QM energies
    diff = sum((mm - qm) ** 2 for mm, qm in zip(relative_mm_energies, qm_energies))
    print(f"Current params: {params}, Difference: {diff}")

    # Number of energy terms from Conformer4-17
    num_terms = 486

    # Calculate the mean absolute error between MM and QM energies
    MAE = sum(np.abs(mm - qm) for mm, qm in zip(relative_mm_energies, qm_energies)) / num_terms
    print(f"MAE: {MAE}")

    # Save the difference to a diff.txt file
    with open('MM_QM_Total_diff_square.txt', 'a') as f:
        f.write(f"Params: {params}, sum_diff_square: {diff}, MAE:{MAE}\n")

    return diff

# Optimization
try:
    folders = ['Conformer3', 'Conformer4', 'Conformer5', 'Conformer6', 'Conformer7', 'Conformer8', 'Conformer9', 'Conformer10', 'Conformer11', 'Conformer12', 'Conformer13', 'Conformer14', 'Conformer15', 'Conformer16', 'Conformer17', 'Conformer18', 'Conformer19', 'Conformer20']

    # Perform optimization
    result = minimize(objective_function, 
                      x0=[ 0.0102, 0.0941, 0.1078, 0.0614, 0.1463, 0.0832, 0.0000], 
                      args=(folders,), 
                      bounds=[(0.0002, 0.2602), (0.0041, 0.3441), (0.0078, 0.3578), (0.0014, 0.3114), 
                              (0.0063, 0.3963), (0.0032, 0.3332), (-2.0000, 0.5000)], 
                      method='L-BFGS-B', 
                      options={'eps': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01],'gtol': 1e-5, 'maxiter': 99999})


    # Save the optimized parameters
    optimized_params = result.x
    with open('optimized_params.json', 'w') as f:
        json.dump(optimized_params.tolist(), f)

    print("Optimized parameters:", optimized_params)
    print("Minimum absolute difference from target:", result.fun)

    # Clean up .txt files
    for folder in folders:
        os.chdir(folder)
        txt_files = [file for file in os.listdir('.') if file.endswith('.txt')]
        for file in txt_files:
            os.remove(file)
        os.chdir('..')

except RuntimeError as e:
    print(f"Optimization failed: {e}")
