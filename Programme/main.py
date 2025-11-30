import sys
import warnings
# On masque les avertissements pour une sortie propre
warnings.filterwarnings("ignore")

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configuration de l'email (OBLIGATOIRE pour NCBI)
Entrez.email = "nathandissezpro@gmail.com"

def importer_sequences():
    """
    Fonction principale pour récupérer les séquences.
    Version corrigée : Feedback détaillé et sans émojis.
    """
    toutes_les_sequences = []
    
    print("\n--- INITIALISATION DE L'IMPORTATION ---")
    
    # 1. Combien d'entrées ?
    try:
        nb_imports = int(input("Combien de sources voulez-vous ajouter ? : "))
    except ValueError:
        print("[ERREUR] Il faut entrer un nombre entier.")
        return []

    # 2. Boucle principale
    for i in range(nb_imports):
        print(f"\n--- SOURCE N° {i + 1} / {nb_imports} ---")
        print("1. NCBI (Numéro d'Accession)")
        print("2. Fichier FASTA local")
        print("3. Saisie manuelle")
        
        choix = input("Votre choix (1/2/3) : ").strip()

        # --- CAS 1 : NCBI ---
        if choix == "1":
            # Liste aide-mémoire dans le prompt
            msg_aide = "Numéro d'accession (ex: M62320, K03455, U46016, K03454, AJ006022, L20587, M15390, M33262, X52154) : "
            acc_id = input(msg_aide).strip()
            
            print(f"[INFO] Connexion au NCBI pour : '{acc_id}'...")
            try:
                # Récupération
                handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="fasta", retmode="text")
                seqs = list(SeqIO.parse(handle, "fasta"))
                handle.close()

                # Vérification
                if len(seqs) > 0:
                    toutes_les_sequences.extend(seqs)
                    print(f"[SUCCES] Séquence récupérée : {seqs[0].id} | Longueur : {len(seqs[0].seq)}")
                else:
                    print(f"[ATTENTION] Aucune séquence trouvée pour l'ID '{acc_id}'. Vérifiez l'orthographe.")

            except Exception as e:
                print(f"[ERREUR] Échec du téléchargement pour '{acc_id}'.")
                print(f"Détail technique : {e}")

        # --- CAS 2 : FICHIER LOCAL ---
        elif choix == "2":
            chemin = input("Chemin du fichier FASTA : ").strip()
            try:
                if list(SeqIO.parse(chemin, "fasta")): # Test rapide si lisible
                    # On recharge proprement
                    seqs = list(SeqIO.parse(chemin, "fasta")) 
                    toutes_les_sequences.extend(seqs)
                    print(f"[SUCCES] {len(seqs)} séquences chargées depuis le fichier.")
                else:
                    print("[ATTENTION] Le fichier est vide ou n'est pas au format FASTA.")
            except Exception as e:
                print(f"[ERREUR] Impossible de lire le fichier '{chemin}'.")
                print(f"Détail technique : {e}")

        # --- CAS 3 : MANUEL ---
        elif choix == "3":
            seq_str = input("Copiez la séquence ici : ").strip().upper()
            nom = input("Nom de la séquence : ").strip()
            
            if seq_str and nom:
                record = SeqRecord(Seq(seq_str), id=nom, description="Import_manuel")
                toutes_les_sequences.append(record)
                print("[SUCCES] Séquence manuelle ajoutée.")
            else:
                print("[ERREUR] Le nom ou la séquence est vide.")

        else:
            print("[ERREUR] Choix invalide. Cette étape est ignorée.")

    print(f"\n--- BILAN : {len(toutes_les_sequences)} SÉQUENCES EN MÉMOIRE ---")
    
    # Petit récapitulatif pour voir ce qui manque
    if len(toutes_les_sequences) > 0:
        print("Liste des séquences valides :")
        for s in toutes_les_sequences:
            print(f" - {s.id}")
    else:
        print("[ATTENTION] Aucune séquence n'a été chargée.")

    return toutes_les_sequences

if __name__ == "__main__":
    importer_sequences()
    