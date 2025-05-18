import os

# Prompt user to choose the metal name
print("Please choose the metal you want to use:")
print("1. MG")
print("2. CA")

metal_choice = input("Enter the number corresponding to your choice: ")

# Map the user's choice to the metal name
metal_map = {"1": "MG", "2": "CA"}

# Validate user input
if metal_choice not in metal_map:
    print("Invalid choice. Please enter a valid number.")
    exit()

Metal = metal_map[metal_choice]

# Prompt user to choose the total time for analysis
while True:
    try:
        total_time_ns = int(input("Please choose the total time in nanoseconds (0-250) you want to use for analysis: "))
        if total_time_ns < 0 or total_time_ns > 250:
            print("Invalid input. Please enter a value between 0 and 250.")
        else:
            break
    except ValueError:
        print("Invalid input. Please enter a valid number.")

# Calculate the frame number
frame = 5001 - int((total_time_ns * 1000) / 50)

# Prompt user to choose the concentration for analysis
print("Please choose the concentration you want to use:")
print("1. 0.05")
print("2. 0.1")
print("3. 0.25")
print("4. 0.5")
print("5. 1")

concentration_choice = input("Enter the number corresponding to your choice: ")

# Map the user's choice to the concentration
concentration_map = {"1": "0.05", "2": "0.1", "3": "0.25", "4": "0.5", "5": "1"}

# Validate user input
if concentration_choice not in concentration_map:
    print("Invalid choice. Please enter a valid number.")
    exit()

Concentration = concentration_map[concentration_choice]

# Create cpptraj input file
with open("cpptraj_vector.in", "w") as f:
    f.write(f"""\
# Specify parameter and trajectory files
parm ../../../entire_box_thf_box.prmtop
trajin ../../../entire_box_thf_box_nvt_prod2.netcdf {frame} last

strip :THF,DMA

""")

    # Count the number of TFS molecules from the pdb file
    with open("../../../entire_box_thf_box.pdb", "r") as pdb_file:
        lines = pdb_file.readlines()
        # Assuming TFS molecules are represented as "TFS" in the file
        tfs_count = sum(1 for line in lines if "TFS" in line)

        # Add vector commands for TFS molecules to cpptraj_vector.in
        for i in range(1, (tfs_count // 15)+1):
            f.write(f"vector center :{i} out TFS{i}_{Concentration}M.dat\n")

        # Add vector commands for the metal molecules to cpptraj_vector.in
        for i in range((tfs_count // 15)+1, ((tfs_count // 5))//2+1):
            f.write(f"vector center :{i} out {Metal}{(i-tfs_count // 15)}_{Concentration}M.dat\n")

    # Add 'run' command to cpptraj_vector.in
    f.write("run\n")

# Run cpptraj with the generated input file
os.system("cpptraj -i cpptraj_vector.in > cpptraj_vector.out")

# Loop through all .dat files in the current directory
for filename in os.listdir("."):
    if filename.endswith(".dat"):
        with open(filename, "r") as f:
            lines = f.readlines()

        # Remove the first line
        lines = lines[1:]

        # Remove columns 5, 6, and 7 and adjust spacing
        modified_lines = []
        for line in lines:
            columns = line.split()
            if len(columns) >= 7:
                modified_line = f"{columns[1]}     {columns[2]}     {columns[3]}\n"
                modified_lines.append(modified_line)

        # Write modified lines back to the file
        with open(filename, "w") as f:
            f.writelines(modified_lines)





