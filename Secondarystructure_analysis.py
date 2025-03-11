import os
import urllib.request
import subprocess
import multiprocessing

def fetch_pdb(pdb_id):
    """Fetch PDB file content directly from RCSB."""
    urls = [
        f"https://files.rcsb.org/download/{pdb_id}.cif",
        f"https://files.rcsb.org/download/{pdb_id}.pdb"
    ]
    for url in urls:
        try:
            with urllib.request.urlopen(url) as response:
                return response.read().decode()
        except Exception as e:
            print(f"Failed to fetch {pdb_id}: {e}")
    return None

def run_stride(pdb_content, pdb_id):
    """Run STRIDE using subprocess and input the PDB content via stdin."""
    try:
        process = subprocess.run(["/Users/hiteshkandarpa/Desktop/Acads/Protein folding/Assignment-1/stride"], 
                                 input=pdb_content, capture_output=True, text=True)
        return process.stdout
    except Exception as e:
        print(f"Error running STRIDE on {pdb_id}: {e}")
        return None

def process_pdb(entry):
    """Fetch PDB content, process it with STRIDE, and save output."""
    pdb_id, chains = entry
    print(f"\n📥 Processing structure '{pdb_id}'...")
    pdb_content = fetch_pdb(pdb_id)
    if pdb_content:
        output = run_stride(pdb_content, pdb_id)
        if output:
            with open("secondary_structure.txt", "a") as f:
                f.write(f">{pdb_id}\n")
                f.write(output + "\n")

def main(pdb_ids_file):
    """Process multiple PDB IDs in parallel."""
    with open(pdb_ids_file, "r") as f:
        lines = f.readlines()[1:]
        lines = [line.strip() for line in lines if line.strip()]
    
    pdb_dict = {}
    for entry in lines:
        fields = entry.split()
        if not fields:
            continue
        pdb_id = fields[0][:4].upper()
        chain_id = fields[0][4].upper() if len(fields[0]) >= 5 else fields[1].upper()
        pdb_dict.setdefault(pdb_id, []).append(chain_id)
    
    pdb_entries = [(pdb_id, list(set(chains))) for pdb_id, chains in pdb_dict.items()]
    
    with multiprocessing.Pool(processes=os.cpu_count()) as pool:
        pool.map(process_pdb, pdb_entries)

if __name__ == "__main__":
    pdb_ids_file = "/Users/hiteshkandarpa/Desktop/Acads/Protein folding/Assignment-1/dunbrack files assignment 1/cullpdb_pc25.0_res0.0-2.0_len40-1000_R0.25_Xray+Nmr_d2025_02_27_chains9858"  # Replace with the actual file path
    main(pdb_ids_file)
