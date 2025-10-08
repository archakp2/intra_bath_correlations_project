#!/bin/env python3

import numpy as np
import os
import subprocess
import tempfile
from datetime import datetime

'''for Top and bottom faces only as detectors'''

rundate = datetime.now().strftime("%Y%m%d_%H%M%S")

# Some global variable for easy tuning

OUTPUT_DIR = "~/myjobs/photons/OQS/"

SCRIPT_NAME = "EntropyCompare.py"

    # short queue has run time less than 6 hours
    # for small n, simulation is fast, so submit to short queue
    # in practice, at n=2 10000 photons take about 3 minutes
    # KEKCC s queue has maximum runtime of 360 minutes


def submit_lsf_job( job_name, command, queue='s', num_cores=1, walltime="5:00", memory="4G"):
    """
    Generates a temporary LSF job script, submits it, and deletes the script after submission.
    
    Parameters:
        job_name (str): Name of the job.
        command (str): Command to execute in the job.
        queue (str): LSF queue to submit to.
        num_cores (int): Number of CPU cores required.
        walltime (str): Maximum runtime in HH:MM format.
        memory (str): Memory requirement (e.g., "4G").
    """
    
    
    # if short then 5 min, if long then 40 min
    walltime = "5:00" if queue=='s' else "40:00"

    # Create a temporary file for the bash script
    #
    with tempfile.NamedTemporaryFile( mode='w', delete=False, suffix='.sh') as script_file:

        script_path = script_file.name
        print( script_path )
        # some extra data
        script_file.write( f"""
#!/bin/bash
#BSUB -J {job_name}
#BSUB -q {queue}
#BSUB -n {num_cores}
#BSUB -W {walltime}
#BSUB -R "rusage[mem={memory}]"
#BSUB -o tmp/{job_name}+{rundate}.out
#BSUB -e tmp/{job_name}+{rundate}.err

cd /home/qup/owaise/myjobs/OQS
{command}
""")

# Job name --> AspRatio_noTIRThreshold_noFresnel + G(number of sides)_n(number of points)_nd(number of direction)_r(raio of the geometry)
# queue --> (either) s (or) l {5:00 hrs or 40:00 hrs walltime}
# walltime --> (depends on the queue) 5:00 (or) 40:00
# memory --> the memory in use
# output file  --> tmp/{job_name}.out
# stdout → normal printed messages from your program (saved in tmp/jobname.out)
    
    # Submit the job using bsub
    try:
        subprocess.run( ["bsub", "-L", "/bin/bash"], shell=True, stdin=open( f"{script_path}", 'r') )
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job: {e}")
    
    # Remove the temporary script file
    
    os.remove(script_path)
    print(f"Job {job_name} submitted and script {script_path} deleted.")




"""
Task here is just to get the python files for various geometry and various refractive indices 
to look out for the resolutin of each

R( geometry, refractive index )
"""
if __name__ == "__main__":

    Number_of_Baths = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])

    for number in Number_of_Baths:
        
        counter = 0 
        n = number + 2  # n = baths + 2

        # --- Bath positions (start from 2) ---
        bathPosns = [i + 2 for i in range(number)]

        # --- Parameters ---
        mu = 0.0
        gamma = 0.01
        beta = 30.0

        mus_str = "-mus " + " ".join([str(mu)] * number)
        gammas_str = "-gamma " + " ".join([str(gamma)] * number)
        betas_str = "-betas " + " ".join([str(beta)] * number)
        bathPos_str = "-bathPos " + " ".join(map(str, bathPosns))

        # --- Dynamic Hamiltonian construction ---
        diag_terms = []
        for i in range(1, n + 1):
            if i == 1:
                val = 0.5
            elif i % 2 == 0:
                val = -0.2
            else:
                val = 0.2
            diag_terms.append(f"h{i}{i} {val}")
        
        coupling_terms = "h21 1 h12 1"
        
        hs_str = "-hs " + " ".join(diag_terms) + " " + coupling_terms

        # --- Example g parameters (modify as needed) ---
        gs_str = "-gs g1 0.1 1.01 0.1"

        # --- Output filenames ---
        suffix = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename1 = f"{OUTPUT_DIR}/B{number:d}_gvsCorelations_t{suffix}.csv"
        filename2 = f"{OUTPUT_DIR}/B{number:d}_gvsC10_t{suffix}.csv"

        # --- Final command ---
        script = (
            f"python3 {SCRIPT_NAME} "
            f"-n {n} {hs_str} {gs_str} "
            f"{gammas_str} {betas_str} {mus_str} "
            f"-W -20 20 {bathPos_str} "
            f"--output1 {filename1} --output2 {filename2}"
        )
    
        counter = counter +1 
        job_name = "OQS_{0:d}".format( counter )
                     
        submit_lsf_job( job_name, script, 'l' )
                    
        counter += 1

    exit(0)
