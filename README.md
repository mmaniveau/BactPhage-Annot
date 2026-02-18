# BactPhage-Annot Pipeline

**BactPhage-Annot** est un pipeline Nextflow DSL2 conçu pour l'annotation de génomes bactériens et des phages. Il intègre une suite d'outils de pointe pour détecter les systèmes de défense, les éléments génétiques mobiles (MGE), les prophages, les intégrons et les systèmes de rétro-évolution (DGR).
L'utilisation des outils est modulable, la suite MGEfinder permet de comparer un échantillon par rapport à une référence pour vérifier l'insertion/délétion de MGE entre les deux.

## 🚀 Fonctionnalités principales

* **Standardisation automatique** : Renommage des contigs (C1, C2...) et annotation via Prokka pour une compatibilité parfaite entre les outils.
* **Systèmes de Défense/AntiDéfense** : Identification des systèmes antiviraux via **DefenseFinder** et des contre-attaques virales via **AntiDefenseFinder**.
* **Éléments Génétiques Mobiles (MGE)** : Détection des insertions via **MGEfinder** (analyse de split-reads) et recherche d'intégrons via **IntegronFinder**.
* **Prophage** : Identification de prophages (**PhiSpy**) et classification des réseaux **CRISPR-Cas**.
* **DGR** : Recherche de Reverse Transcriptases et de systèmes **DGR** (Diversity-Generating Retroelements).
* **Synthèse de données** : Génération d'un tableau CSV récapitulatif et d'un fichier GFF fusionné.

---

## 📋 Pré-requis

* **Nextflow** (DSL2)
* **Docker** à télécharger/monter
* **Modèles DefenseFinder** : Doivent être placés dans le dossier `defensefinder_models/` à la racine

## 📄 Structure du fichier d'entrée (`samples.csv`)

Le pipeline nécessite un fichier CSV nommé par défaut `samples.csv`. Ce fichier fait le lien entre vos identifiants d'échantillons et les différents fichiers bruts ou pré-annotés nécessaires à chaque outil.

### 📋 En-têtes et colonnes
Le fichier doit comporter les 5 colonnes suivantes (séparées par des virgules) :

| Colonne | Description | Usage principal |
| :--- | :--- | :--- |
| `sample_id` | Identifiant unique de l'échantillon | Nommage des dossiers/fichiers |
| `fasta` | Assemblage génomique (contigs) | CRISPR, DGR, Integrons, Defense |
| `bam` | Alignement des reads de l'échantillon | MGEfinder (détection split-reads) |
| `bai` | Index du fichier BAM | MGEfinder |
| `ref_fasta` | Génome de référence (ex: souche ancestrale) | MGEfinder (comparaison InDel) |

### 💡 Exemple de contenu
Voici à quoi doit ressembler votre fichier (avec des chemins absolus recommandés) :

``csv
sample_id,fasta,bam,bai,ref_fasta
Sp_1,/Sp_1.fasta,/Sp_1.bam,/Sp_1.bam.bai,/Ancestral.fasta
Sp_2,/Sp_2.fasta,/Sp_2.bam,/Sp_2.bam.bai,/Ancestral.fasta

## ⚙️ Paramètres de Configuration (nextflow.config)

Le pipeline est configurable via le bloc `params`. Voici le détail des options disponibles :

### 🛠️ Activation des modules
| Paramètre | Description | Défaut |
| :--- | :--- | :--- |
| `run_mgefinder` | Active la détection d'éléments mobiles (BAM requis) | `true` |
| `run_integron` | Active la recherche d'intégrons | `true` |
| `run_phispy` | Active la détection de prophages (GBK requis) | `true` |
| `run_defense` | Active DefenseFinder et AntiDefenseFinder | `true` |
| `run_crispr` | Active CRISPRCasFinder | `true` |
| `run_dgr` | Active la recherche de DGR | `true` |

### 🧬 Paramètres spécifiques aux outils

#### **IntegronFinder**
* `min_attc_size` / `max_attc_size` : Plage de taille (bp) pour les sites attC (Défaut: `40-200`).
* `distance_thresh` : Distance max (bp) pour lier des éléments entre eux (`4000`).
* `topology` : Structure de l'ADN (`linear` ou `circ`).

#### **PhiSpy**
* `phispy_threads` : Nombre de cœurs CPU alloués (`2`).
* `phispy_min_contig_size` : Taille minimale des contigs à analyser (`10`).
* `phispy_window_size` : Fenêtre glissante pour l'algorithme (`30`).

#### **DefenseFinder**
* `defense_coverage` : Couverture minimale des profils HMM (`0.4`).
* `defense_db_type` : Mode d'organisation génomique (`ordered_replicon`).
* `defensefinder_models` : Chemin local vers les modèles HMM.

#### **CRISPRCasFinder**
* `crispr_vicinity` : Fenêtre de recherche (bp) autour des arrays pour les gènes Cas (`20000`).
* `crispr_cas` : Active la recherche spécifique des gènes Cas (`true`).


## 🐳 Gestion des Conteneurs (Docker)

Le pipeline utilise **Docker** pour garantir la reproductibilité. Les images sont configurées automatiquement pour chaque processus :

* **MGEfinder** : `ivasilyev/mgefinder` | DockerHub |
* **PhiSpy** : `drantoine/phispy:1.0.2` | DockerHub |
* **IntegronFinder** : `gempasteur/integron_finder:2.0.6` | DockerHub |
* **DefenseFinder** : `bjhall/defense-finder:2.0.1-models-2.0.2` | DockerHub |
* **AntiDefenseFinder** : `bjhall/defense-finder:2.0.1-models-2.0.2` | DockerHub |
* **CRISPRCasFinder** : `unlhcc/crisprcasfinder:4.2.20` | DockerHub |
* **DGR** : `dgr:1.0` | **Custom Dockerfile** | docker build -t dgr:1.0 -f dgr1.0.dockerfile .

Comme il n'existe pas d'image officielle pour **DGR-package**, un conteneur dédié a été conçu (basé sur Ubuntu 22.04). Il assure la présence de :
* **DGR-package** (via pip) : Le cœur du module pour l'identification des Reverse Transcriptases et l'alignement TR/VR.
* **Outils Bioinformatiques** : `prodigal` (prédiction d'ORFs), `hmmer` (recherche de domaines RT), et `ncbi-blast+` (alignements locaux).
* **Dépendances Python** : `biopython`, `pandas` et `DGR-package`.

> **Note importante** : Pour DefenseFinder, le pipeline monte automatiquement un volume local (`-v`) pour accéder aux modèles de défense situés dans votre répertoire de projet. Le Docker utilise des modèles Crispr 2.0 : https://github.com/macsy-models/CasFinder/releases/tag/3.1.0

# 🖥️ Utilisation du Pipeline

Pour lancer l'analyse, utilisez la commande suivante. L'activation de Docker est obligatoire pour garantir le bon fonctionnement des outils.

### Exemple de commande brute
```bash
nextflow run main.nf --sample_csv /user/home/pipeline/CSV/samples.csv -with-docker --run_name run_01 --run_mgefinder false
```
--sample_csv : Chemin absolu vers votre fichier CSV d'entrée.

-with-docker : Active les conteneurs (requis).

--run_name : Nom du dossier de sortie.

**ATTENTION** : Si vous relancez cette commande sans changer le nom, le dossier existant sera écrasé.

--run_mgefinder false : Désactive le module MGEfinder (utile si vous n'avez pas de données BAM/référence).

Si le pipeline est interrompu, ajoutez -resume pour ne pas recalculer ce qui a déjà été fait.

## 📂 Organisation des Résultats

Après l'exécution, les résultats sont organisés comme suit :

* **`_Final_Summary/`** :
    * `BactPhage-Annot_summary.csv` : Tableau croisé de tous les résultats par échantillon.
* **`_GFF_Summary/`** :
    * `*_MGE_map.gff` : Fichier GFF3 fusionné contenant toutes les annotations détectées (visualisable dans IGV ou Artemis).
* **Dossiers par outil** :
    * `defense/`, `anti_defense/`, `phispy/`, `integron/`, `mge_find/`, `crispr/`, `dgr/`.

---

## 🔬 Détails des Outils

### MGE & Intégrons
* **MGEfinder** : Identifie les insertions en comparant un fichier d'alignement (**BAM**) à une référence.

La suite **MGEfinder** permet une analyse comparative entre un échantillon et une référence afin d'identifier les événements d'**insertion** ou de **délétion**. 

> **⚠️ Note importante sur la préparation des fichiers BAM :**
> Pour que MGEfinder fonctionne correctement, vos fichiers d'entrée doivent respecter les critères suivants :
> * **Indexation** : Chaque fichier `.bam` doit être accompagné de son index `.bai` (ex: `mon_fichier.bam` et `mon_fichier.bam.bai`).
> * **Tri** : Les fichiers BAM doivent impérativement être **triés par coordonnées** (`samtools sort`).
> * **Tags MD** : Les fichiers doivent contenir le tag `MD` (généré par `bwa mem`). Si vos fichiers BAM ont été produits par un autre aligneur, vous devez les traiter au préalable avec la commande suivante :
>   `samtools calmd -b mon_fichier.bam ref_fasta > mon_fichier_avec_md.bam`

### Intégrons
* **IntegronFinder** : Recherche les sites de recombinaison **attC** et les gènes d'**intégrases**.

### Prophage
* **PhiSpy** : Détecte les régions prophagiques en combinant des méthodes statistiques (**Random Forest**) et des données de similarité. Il analyse des critères comme le taux de GC, la taille des gènes et les blocs de transcription unidirectionnels pour identifier de nouveaux phages.

### Défense Antivirale
* **DefenseFinder** : Utilise des profils **HMM** et **MacSyFinder** pour identifier et valider les systèmes complets par analyse de synténie.


### CRISPR & DGR
* **CRISPRCasFinder** : Détecte les *arrays* (répétitions/spacers) et les gènes **Cas** associés.

* **DGRpy** : Identifie les **Reverse Transcriptases** et valide les alignements **TR/VR** .

---

## ⚠️ Gestion des exceptions

Le pipeline inclut une gestion robuste des erreurs spécifiques pour garantir une exécution sans interruption :

* **PhiSpy** : Les codes d'erreur **30** (contigs trop courts) et **41** (manque d'ORFs) sont capturés. En cas d'échec, le pipeline génère un fichier `README_phispy_skipped.txt` plutôt que de s'arrêter.
* **Fichiers vides** : Les scripts de synthèse ignorent automatiquement les fichiers corrompus ou vides, assurant ainsi la production du rapport final même en cas de résultats nuls pour certains outils.

## Citations

* **MGEfinder**

Durrant, M. G., Li, M. M., Siranosian, B. A., Montgomery, S. B. & Bhatt, A. S. A Bioinformatic Analysis of Integrative Mobile Genetic Elements Highlights Their Role in Bacterial Adaptation. Cell Host & Microbe 0, (2019)

* **IntegronFinder**

The paper is published in Microorganism.

Néron, Bertrand, Eloi Littner, Matthieu Haudiquet, Amandine Perrin, Jean Cury, and Eduardo P.C. Rocha. 2022. IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella Microorganisms 10, no. 4: 700. https://doi.org/10.3390/microorganisms10040700

Nawrocki, E.P. and Eddy, S.R. (2013) Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29, 2933-2935.
Eddy, S.R. (2011) Accelerated Profile HMM Searches. PLoS Comput Biol, 7, e1002195.
Hyatt, D., Chen, G.L., Locascio, P.F., Land, M.L., Larimer, F.W. and Hauser, L.J. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11, 119.

and if you use the function --func_annot which uses NCBIfam-AMRFinder hmm profiles:

    Haft, DH et al., Nucleic Acids Res. 2018 Jan 4;46(D1):D851-D860 PMID: 29112715


* **PhiSpy**
  
Akhter et al

Sajia Akhter, Ramy K. Aziz, Robert A. Edwards; PhiSpy: a novel algorithm for finding prophages in bacterial genomes that combines similarity- and composition-based strategies. Nucl Acids Res 2012; 40 (16): e126. doi: 10.1093/nar/gks406

For the newest versions, please also cite:

McNair, K., Decewicz, P., Akhter, S., Aziz, R.K., Daniel, S., Edwards, R.A. 2019. PhiSpy. https://github.com/linsalrob/PhiSpy/ doi://10.5281/zenodo.3475717

* **DefenseFinder/AntiDefenseFinder**
  
"Systematic and quantitative view of the antiviral arsenal of prokaryotes" Nature Communication, 2022, Tesson F., Hervé A. , Mordret E., Touchon M., d’Humières C., Cury J., Bernheim A.

"MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes." Peer Community Journal, Volume 3 (2023), article no. e28. Néron, Bertrand; Denise, Rémi; Coluzzi, Charles; Touchon, Marie; Rocha, Eduardo P.C.; Abby, Sophie S.

"CRISPRCasFinder, an update of CRISRFinder, includes a portable version, enhanced performance and integrates search for Cas proteins." Nucleic Acids Research 2018 Couvin D. et al.

* **CRISPRCasFinder**

Grissa I, Vergnaud G, Pourcel C. CRISPRFinder: a web tool to identify clustered regularly interspaced short palindromic repeats. Nucleic Acids Res. 2007 Jul;35(Web Server issue):W52-7. DOI: https://doi.org/10.1093/nar/gkm360 PMID:17537822

Abby SS, Néron B, Ménager H, Touchon M, Rocha EP. MacSyFinder: a program to mine genomes for molecular systems with an application to CRISPR-Cas systems. PLoS One. 2014 Oct 17;9(10):e110726. DOI: https://doi.org/10.1371/journal.pone.0110726 PMID:25330359

Couvin D, Bernheim A, Toffano-Nioche C, Touchon M, Michalik J, Néron B, Rocha EPC, Vergnaud G, Gautheret D, Pourcel C. CRISPRCasFinder, an update of CRISRFinder, includes a portable version, enhanced performance and integrates search for Cas proteins. Nucleic Acids Res. 2018 Jul 2;46(W1):W246-W251. DOI: https://doi.org/10.1093/nar/gky425 PMID:29790974

Néron B, Denise R, Coluzzi C, Touchon M, Rocha EPC, Abby SS (2022). MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes. bioRxiv DOI: 10.1101/2022.09.02.506364 

Further information are available at: https://crisprcas.i2bc.paris-saclay.fr.

* **DGR**

https://pypi.org/project/DGR-package/

Roux, S., Paul, B.G., Bagby, S.C. et al. Ecology and molecular targets of hypermutation in the global microbiome. Nat Commun 12, 3076 (2021). https://doi.org/10.1038/s41467-021-23402-7

Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119

Eddy SR (2011) Accelerated Profile HMM Searches. PLoS Comput Biol 7(10): e1002195. https://doi.org/10.1371/journal.pcbi.1002195

* **Prokka**
  
Torsten Seemann, Prokka: rapid prokaryotic genome annotation, Bioinformatics, Volume 30, Issue 14, July 2014, Pages 2068–2069, https://doi.org/10.1093/bioinformatics/btu153
