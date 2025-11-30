import sys
import os
import warnings
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# On masque les avertissements pour la soutenance
warnings.filterwarnings("ignore")

# Configuration Email (Obligatoire pour NCBI)
Entrez.email = "nathandissezpro@gmail.com"

def extraire_genes_interessants(record_genbank):
    """
    Fonction : Reçoit une fiche GenBank complète et extrait Env et Tat.
    Retourne un dictionnaire avec les séquences découpées.
    """
    genes_trouves = {"env_prot": None, "tat_prot": None}
    
    # On parcourt toutes les caractéristiques/features
    for feature in record_genbank.features:
        if feature.type == "CDS": # CDS = Coding Sequence
            
            # On récupère le nom du gène
            nom_gene = feature.qualifiers.get("gene", [""])[0].lower()
            descript_produit = feature.qualifiers.get("product", [""])[0].lower()
            
            # --- ENV ---
            if "env" in nom_gene or "envelope" in descript_produit:
                # On prend la traduction en protéine déjà faite par NCBI
                if "translation" in feature.qualifiers:
                    genes_trouves["env_prot"] = feature.qualifiers["translation"][0]
            
            # --- TAT ---
            elif "tat" in nom_gene:
                if "translation" in feature.qualifiers:
                    genes_trouves["tat_prot"] = feature.qualifiers["translation"][0]

    return genes_trouves

def pipeline_principale():
    print("\n=== PIPELINE D'ACQUISITION ET EXTRACTION ===")
    
    # virus cibles de la consigne
    cibles = ["M62320", "K03455", "U46016", "K03454", "AJ006022", "L20587", "M15390", "M33262", "X52154"]
    
    sequences_env = []
    sequences_tat = []
    
    print(f"[INFO] Demarrage du telechargement pour {len(cibles)} souches virales...")

    for acc_id in cibles:
        try:
            # 1. TÉLÉCHARGEMENT (Format GenBank pour avoir les coordonnées)
            handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # 2. EXTRACTION
            resultats = extraire_genes_interessants(record)
            
            # 3. STOCKAGE
            nom_virus = f"{acc_id}" # On garde l'ID comme identifiant
            
            # Si on trouve Env
            if resultats["env_prot"]:
                seq_obj = SeqRecord(Seq(resultats["env_prot"]), id=nom_virus, description="Env_Protein")
                sequences_env.append(seq_obj)
            else:
                print(f"[ATTENTION] Pas de gene Env trouve pour {acc_id}")

            # Si on trouve Tat
            if resultats["tat_prot"]:
                seq_obj = SeqRecord(Seq(resultats["tat_prot"]), id=nom_virus, description="Tat_Protein")
                sequences_tat.append(seq_obj)
            else:
                print(f"[ATTENTION] Pas de gene Tat trouve pour {acc_id}")

            print(f"[SUCCES] Traite : {acc_id}")
            
        except Exception as e:
            print(f"[ERREUR] Probleme sur {acc_id} : {e}")

    # 4. SAUVEGARDE DES FICHIERS DANS LE SOUS-DOSSIER
    print("\n--- GENERATION DES FICHIERS ---")
    
    # Definition du dossier de stockage (Un etage plus haut, puis 'Donnees')
    dossier_cible = "../Donnees"
    
    # On cree le dossier s'il n'existe pas deja (exist_ok=True evite de planter si deja la)
    if not os.path.exists(dossier_cible):
        os.makedirs(dossier_cible)
        print(f"[INFO] Dossier '{dossier_cible}' cree.")

    if sequences_env:
        # On joint le dossier et le nom de fichier
        chemin_env = os.path.join(dossier_cible, "env_all.fasta")
        SeqIO.write(sequences_env, chemin_env, "fasta")
        print(f"[SUCCES] Fichier genere : '{chemin_env}' ({len(sequences_env)} sequences)")
    
    if sequences_tat:
        chemin_tat = os.path.join(dossier_cible, "tat_all.fasta")
        SeqIO.write(sequences_tat, chemin_tat, "fasta")
        print(f"[SUCCES] Fichier genere : '{chemin_tat}' ({len(sequences_tat)} sequences)")
        
if __name__ == "__main__":
    pipeline_principale()