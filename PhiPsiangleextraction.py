import os
import urllib.request
import pickle  # ✅ Import pickle for saving data
from Bio import PDB
from io import StringIO
import math
from Bio.PDB.Polypeptide import is_aa, PPBuilder, Polypeptide
import multiprocessing as mp

def fetch_structure(pdb_id):
    i=0
    urls = [
        f"https://files.rcsb.org/download/{pdb_id}.cif",
        f"https://files.rcsb.org/download/{pdb_id}.pdb"
    ]
    for url in urls:
        try:
            response = urllib.request.urlopen(url)
            pdb_data = response.read().decode()
            i+=1
            print(f"✅ Successfully fetched: {url} Number - {i}, {9859-i} to go")
            return StringIO(pdb_data)
        except Exception as e:
            print(f"⚠️ Failed to fetch: {url} ({e})")
    print(f"❌ Error: Structure {pdb_id} not found in either CIF or PDB format.")
    return None

def get_phi_psi_angles(pdb_id, structure_data):
    parser = PDB.MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, structure_data)
        model = structure[0]
        ppb = PPBuilder()
        all_angles = {}

        for chain in model:
            chain_id = chain.id
            peptides = ppb.build_peptides(chain)
            if not peptides:
                print(f"⚠️ No peptides found in Chain {chain_id} of {pdb_id}. Trying fallback...")
                residues = [res for res in chain if is_aa(res, standard=True)]
                if residues:
                    poly = Polypeptide(residues)
                    peptides = [poly]
                else:
                    print(f"❌ No valid amino acid residues in {pdb_id} Chain {chain_id}.")
                    continue

            angles_list = []
            for pp in peptides:
                phi_psi_angles = pp.get_phi_psi_list()
                if not phi_psi_angles:
                    print(f"❌ No phi/psi angles found for {pdb_id} Chain {chain_id}.")

                for residue, (phi, psi) in zip(pp, phi_psi_angles):
                    res_name = residue.get_resname()
                    res_number = residue.get_id()[1]
                    phi_val = round(math.degrees(phi), 2) if phi is not None else "N/A"
                    psi_val = round(math.degrees(psi), 2) if psi is not None else "N/A"

                    angles_list.append((res_name, res_number, phi_val, psi_val))
                    print(f"📌 {pdb_id} Chain {chain_id} Residue {res_number} ({res_name}): Phi={phi_val}, Psi={psi_val}")

            if angles_list:
                all_angles[chain_id] = angles_list

        return all_angles
    except Exception as e:
        print(f"❌ Error parsing structure {pdb_id}: {e}")
        return {}

def process_pdb(entry):
    pdb_id, desired_chains = entry
    print(f"\n📥 Processing structure '{pdb_id}'...")
    structure_data = fetch_structure(pdb_id)
    if not structure_data:
        return (pdb_id, {})
    all_angles = get_phi_psi_angles(pdb_id, structure_data)
    if desired_chains:
        filtered = {chain: angles for chain, angles in all_angles.items() if chain in desired_chains}
    else:
        filtered = all_angles
    return (pdb_id, filtered)

if __name__ == '__main__':
    file_path = "/Users/hiteshkandarpa/Desktop/Acads/Protein folding/dunbrack files assignment 1/cullpdb_pc25.0_res0.0-2.0_len40-1000_R0.25_Xray+Nmr_d2025_02_27_chains9858"  # Change to your file path

    with open(file_path, "r") as f:
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

    pool = mp.Pool(mp.cpu_count())
    results_list = pool.map(process_pdb, pdb_entries)
    pool.close()
    pool.join()

    final_results = {pdb_id: result for pdb_id, result in results_list}

    # ✅ Save final_results to a file
    with open("final_results.pkl", "wb") as f:
        pickle.dump(final_results, f)

    print("✅ Saved final_results to 'final_results.pkl'")

    # ✅ Print summary
    print("\n📌 **Extracted Phi/Psi Angles**\n")
    header = f"{'PDB ID':<8} {'Chain':<6} {'Residue':<8} {'AA':<4} {'Phi (°)':<10} {'Psi (°)':<10}"
    print(header)
    print("=" * len(header))

    for pdb_id, chains in final_results.items():
        for chain_id, angles in chains.items():
            if not angles:
                print(f"{pdb_id:<8} {chain_id:<6} ❌ No phi/psi angles found.")
            else:
                for res_name, res_number, phi, psi in angles:
                    print(f"{pdb_id:<8} {chain_id:<6} {res_number:<8} {res_name:<4} {phi:<10} {psi:<10}")

    print("\n✅ **Extraction complete!**")
