#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -------------------------------
// Lecture et validation CSV
// -------------------------------

//Vérification de l'existence du CSV
def csvFile = file(params.sample_csv)
if (!csvFile.exists()) {
    error "[ERROR] Fichier CSV '${params.sample_csv}' introuvable."
}
println "[INFO] CSV fourni: ${params.sample_csv}"
println "[INFO] Output directory: ${params.outdir}\n"

/*
Création des channels pour les process
Le pipeline lit le CSV ligne par ligne en vérifiant l'existence de chaque paramètre
*/
samples_ch = Channel
    .fromPath(params.sample_csv)
    .splitCsv(header: true)
    .map { row ->

        if (!row.sample_id || !row.fasta) {
            error "[ERROR] Colonnes obligatoires manquantes (sample_id, fasta)"
        }

        def sample_id = row.sample_id.trim()

        // Vérification Fasta
        def fasta = file(row.fasta)
        if (!fasta.exists()) {
            error "[ERROR] FASTA introuvable pour '${sample_id}' : ${row.fasta}"
        }

        // Vérification BAM / BAI / REF (MGEfinder)
        def bam       = row.bam && row.bam.trim() ? file(row.bam) : null
        def bai       = row.bai && row.bai.trim() ? file(row.bai) : null
        def ref_fasta = row.ref_fasta && row.ref_fasta.trim() ? file(row.ref_fasta) : null

        if (params.run_mgefinder) {
            if (!bam || !bam.exists()) error "[ERROR] BAM introuvable pour '${sample_id}' : ${row.bam}"
            if (!bai || !bai.exists()) error "[ERROR] BAI introuvable pour '${sample_id}' : ${row.bai}"
            if (!ref_fasta || !ref_fasta.exists()) error "[ERROR] FASTA introuvable pour '${sample_id}' : ${row.ref_fasta}"
        }

        // Retourne un tuple contenant l'ensemble des paramètres
        return [sample_id, fasta, bam, bai, ref_fasta]
    }

process PREPARE_DATA {
    tag "${sample_id}"
    container 'staphb/prokka:1.14.6'

    input:
    tuple val(sample_id), path(fasta), path(bam), path(bai), path(ref)

    output:
    tuple val(sample_id), path("clean.fasta"), path("clean.gff"), path(bam), path(bai), path(ref), path("clean.gbk"), path("${sample_id}_mapping.tsv"), emit: ready

    script:
    """
    # 1. Création de la table de correspondance à 3 colonnes
    # On extrait les headers, et pour chaque ligne on affiche :
    # Nom_Court <tab> Nom_Original <tab> Sample_ID
    grep "^>" ${fasta} | sed 's/^>//' > names.txt
    awk -v sam="${sample_id}" 'BEGIN{FS=OFS="\\t"} {print "C"NR, \$0, sam}' names.txt > ${sample_id}_mapping.tsv

    # 2. Renommage du FASTA pour la compatibilité (C1, C2...)
    awk '/^>/ {print ">C" ++i; next} {print}' ${fasta} > renamed.fasta

    # 3. Prokka pour générer un GenBank propre
    prokka --outdir out_prokka --prefix clean --locustag C --cpus 8 --mincontiglen 0 renamed.fasta
    
    mv out_prokka/clean.fna ./clean.fasta
    mv out_prokka/clean.gbk ./clean.gbk
    cp out_prokka/clean.gff ./clean.gff
    """
}

// -------------------------------
// Processus : MGEfinder find
// -------------------------------

/*
DESCRIPTION
 MGEfinder est un pipeline de détection de MGE (élément génétique mobile) séparé en 3 modules. 
 Le premier "fin" scanne le fichier d'alignement (BAM) pour identifier des "discordant reads" (reads mal alignés ou split-reads) 
 qui signalent la présence potentielle d'un MGE. 
 Si une bactérie par exemple a acquis un nouvel élément mobile que la référence n'a pas,
 les reads peuvent être coupés (split-reads) 

GÉNOME (---) | MGE (___)
             |
-------|-----|----------    <-- Ton fragment d'ADN (150 bp)
[ R1 >>>>>>>>>> ]           <-- Read 1 (100 bp) : 50bp génome + 50bp MGE
        [ <<<<<<<<<< R2 ]   <-- Read 2 (100 bp) : 50bp MGE + 50bp génome
             |
          JONCTION


REF : ---------------------------------
             |
R1 :  [ >>>>>|XXXXX]        <-- "Split-read" !
             |
R2 :         [XXXXX|<<<<< ] <-- "Split-read" aussi !
             |
             |--> Position X : MGEfinder note "CANDIDAT" ici.


INPUT
- sample_id : Identifiant de l'échantillon
- ref_fasta : Le génome de référence qui servira pour la suite de MGEfinder
- bam       : Fichier d'alignement trié par coordonnées
            --> Les fichiers BAM d'entrée doivent être triés par coordonnées et indexés
            --> Ils doivent contenir les tags MD (générés par BWA MEM)
                ou traités avec samtools calmd' pour ajouter les tags MD manquants
- bai       : Index du fichier BAM fourni par l'utilisateur
            --> Format : mon_fichier.bam.bai 
                selon le fichier bam correspondant : mon_fichier.bam

OUTPUT 
- Un fichier TSV (${sample_id}.find.tsv) listant tous les sites d'insertions candidats
- Le tuple d'entrée est transféré au module MGEfinder "pair" pour maintenir une cohérence entre module
*/

process mge_find {

    tag "$sample_id"

    publishDir "${params.outdir}/${params.run_name}/mge_find/${sample_id}", mode: 'copy'

    container 'ivasilyev/mgefinder@sha256:4a9942869fe17d2543b73af75ca5ed96f32e35e9de51ef011faeb123ea619efa'

    input:
    tuple val(sample_id), path(ref_fasta), path(bam), path(bai)

    output:
    tuple val(sample_id), path(ref_fasta), path(bam), path("${sample_id}.find.tsv"), path(bai)

    script:
    """
    # Indexer la référence pour BWA
    bwa index ${ref_fasta}

    mgefinder find -id ${sample_id} -o ${sample_id}.find.tsv ${bam}

    """
}


// -------------------------------
// Processus : MGEfinder pair
// -------------------------------

/*
DESCRIPTION
 MGEfinder est un pipeline de détection de MGE (élément génétique mobile) séparé en 3 modules. 
 Le second module "pair" analyse les candidats détectés par "find" pour regrouper les signaux 
 provenant des deux extrémités (5' et 3') d'un même élément mobile. 
 Il "apparie" les jonctions pour définir les bornes précises de l'insertion.

INPUT
- sample_id     : Identifiant de l'échantillon
- ref_fasta     : Le génome de référence
- bam           : Fichier d'alignement (utilisé pour extraire les séquences des jonctions)
- find_tsv      : Résultat du module MGEfinder "find" contenant les sites candidats
- bai           : Index du fichier BAM fourni par l'utilisateur

OUTPUT 
- Un fichier TSV (${sample_id}.pair.tsv) contenant les paires de jonctions validées.
- Le tuple de sortie est transféré au module final "infer" pour la reconstruction des séquences.
*/

process mge_pair {

    tag "$sample_id"

    publishDir "${params.outdir}/${params.run_name}/mge_pair/${sample_id}", mode: 'copy'

    container 'ivasilyev/mgefinder@sha256:4a9942869fe17d2543b73af75ca5ed96f32e35e9de51ef011faeb123ea619efa'

    input:
    tuple val(sample_id), path(ref_fasta), path(bam), path(find_tsv), path(bai)

    output:
    tuple val(sample_id), path(ref_fasta), path("${sample_id}.pair.tsv")

    script:
    """
    # On vérifie si find_tsv contient des candidats
    if grep -qv "^sample_id" ${find_tsv} ; then
        mgefinder pair ${find_tsv} ${bam} ${ref_fasta} > ${sample_id}.pair.tsv
    else
        # On crée un fichier pair vide avec l'en-tête attendu par inferseq sinon il contiendra des logs
        echo -e "sample_id\tpair_id\tseq_5p\tseq_3p\tcontig\tsite\tjunction_type" > ${sample_id}.pair.tsv
    fi
    """
}


// -------------------------------------
// Processus : MGEfinder inferseq
// -------------------------------------

/*
DESCRIPTION
 MGEfinder est un pipeline de détection de MGE (élément génétique mobile) séparé en 3 modules. 
 Le troisième module "inferseq" est l'étape de reconstruction. Il utilise les paires de 
 jonctions validées par le module "pair" et le génome de référence pour déduire la 
 séquence nucléotidique complète de l'élément inséré.

INPUT
- sample_id : Identifiant de l'échantillon
- ref_fasta : Le génome de référence
- pair_tsv  : Résultat du module MGEfinder "pair" contenant les jonctions appariées

OUTPUT 
- Un fichier TSV (${sample_id}.inferseq.tsv) contenant les séquences reconstruites des MGE
  --> Ce fichier est le résultat final exploitable pour l'annotation ou l'analyse comparative
*/

process mge_inferseq {

    tag "$sample_id"

    publishDir "${params.outdir}/${params.run_name}/mge_inferseq/${sample_id}", mode: 'copy'

    container 'ivasilyev/mgefinder@sha256:4a9942869fe17d2543b73af75ca5ed96f32e35e9de51ef011faeb123ea619efa'

    input:
    tuple val(sample_id), path(ref_fasta), path(pair_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.inferseq.tsv"), emit: infer_out

    script:
    """
    # On regarde la première ligne. 
    # Si elle commence par 'sample_id', c'est un vrai TSV.
    # Si elle commence par '#', ce sont des logs de mgefinder qui n'a rien trouvé.
    if grep -q "^sample_id" ${pair_tsv} ; then
        mgefinder inferseq-reference ${pair_tsv} ${ref_fasta} > ${sample_id}.inferseq.tsv
    else
        # On crée un fichier de sortie vide compatible avec la suite
        echo -e "sample_id\tpair_id\tseq" > ${sample_id}.inferseq.tsv
    fi
    """
}

// -------------------------------
// IntegronFinder
// -------------------------------

/*
DESCRIPTION
IntegronFinder recherche les intégrons mobiles et chromosomiques. 
Il détecte les gènes essentiels (integrase intI, int2..), 
les sites de recombinaison (attC) et les cassettes de gènes associées.

PARAMÈTRES RÉGLABLES (via nextflow.config)
- params.min_attc_size / max_attc_size : Taille (bp) attendue des sites attC (Défaut: 40-200)
- params.distance_thresh       : Distance max entre deux éléments pour les lier (Défaut: 4000bp)
- params.calin_threshold     : Nombre minimum de sites attC pour former un CALIN (Défaut: 2)
- params.evalue_attc         : Seuil statistique pour la détection des sites attC (Défaut: 1)
- params.topology            : Structure de l'ADN ('linear' ou 'circ') pour la détection aux jonctions

INPUT
- sample_id : Identifiant de l'échantillon
- fasta     : Séquences nucléotidiques des contigs

OUTPUT 
- Un dossier (${sample_id}/) contenant tous les résultats de l'outil :
  --> .integrons : Fichier principal listant les intégrons détectés (complets, incal ou CALIN)
    --> complet : + une intégrase + un site attC minimum 
    + l'intégrase et le premier attC doivent être dans une fenêtre proche
    pour être considérer comme liés fonctionnellement (--distance-thresh, par défaut 4000bp)
    --> incal : + une intégrase - sans site attC
    --> calin : - pas d'intégrase + un ou plusieurs site attC regroupés
  --> .summary   : Résumé statistique des éléments trouvés
  --> Fichiers GenBank et fichiers de séquences des cassettes
*/

process run_integronFinder {
    
    tag "${sample_id}"

    publishDir "${params.outdir}/${params.run_name}/integron", mode: 'copy'

    container 'gempasteur/integron_finder:2.0.6'
    // Indispensable car le conteneur ne peut utiliser que les commandes integronFinder
    // Sans ça, il ne reconnaît pas les commandes Bash (mkdir, mv) 
    // et renvoie une erreur "-c eval".
    containerOptions '--entrypoint ""'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: int_results

    script:
    """
    integron_finder --local-max \
        --func-annot \
        --gbk \
        --min-attc-size ${params.min_attc_size} \
        --max-attc-size ${params.max_attc_size} \
        --distance-thresh ${params.distance_thresh} \
        --calin-threshold ${params.calin_threshold} \
        --evalue-attc ${params.evalue_attc} \
        --${params.topology} \
        --outdir out_integron \
        ${fasta}

    #Par défaut integronFinder met tous les résultats d'échantillon dans le même dossier
    #Ici les résultats sont séparés dans un dossier par échantillon
    mkdir -p ${sample_id}

    #Simplifie le nom des dossiers en utilisant seulement le sample_id
    mv out_integron/Results_Integron_Finder_*/* ${sample_id}/
    """
}

// -------------------------------
// PhiSpy
// -------------------------------

/*
DESCRIPTION
PhiSpy est un outil de détection de prophages.
Il combine des méthodes statistiques (random forest) et des
données de similarité de gènes (alignement) pour identifier les régions virales.
Alorithme Random Forest selon différents critères 
(taille des gènes, bloc de transcription unidirectionnel, 
taux de GC, détection de phage like/hypothetical protein,k-mer...)
--> Justifie l'utilisation d'un .gbk pour avoir une annotation des gènes
--> Machine learning qui recherche les signaux correspondant 
à un prophage permettant la détection de potentiels nouveaux
prophages
--> Combinaison avec un programme d'alignement pour comparer
à des virus connus

PARAMÈTRES RÉGLABLES (via nextflow.config)
- params.phispy_threads         : Nombre de cœurs CPU alloués (Défaut: 1)
- params.phispy_min_contig_size : Taille minimale (bp) pour analyser un contig (Défaut: 2000)
- params.phispy_window_size     : Taille de la fenêtre glissante pour l'analyse (Défaut: 30)
- params.phispy_number          : Nombre minimum d'ORFs viraux pour valider un prophage (Défaut: 5)
- params.phispy_extra_dna       : Longueur de l'ADN à inclure autour du prophage (Défaut: 2000)
- params.phispy_min_repeat_len  : Longueur min des répétitions aux extrémités (Défaut: 10)

INPUT
- sample_id : Identifiant de l'échantillon
- gbk       : Fichier GenBank complet (nécessaire car PhiSpy utilise les annotations 
              pour calculer les ratios de gènes viraux/bactériens)

OUTPUT 
- Un dossier (${sample_id}/) contenant tous les résultats de l'outil :
  --> phispy.log                : Log d'exécution de la commande par échantillon
  --> prophage_coordinates.tsv  : Fichier principal listant les coordonnées des prophages
  --> prophage.gff3             : Annotation des prophages pour visualisation (IGV par exemple)
  --> README_phispy_skipped     : Fichier créé manuellement si l'échantillon est sauté 
                                (car contigs trop courts ou manque d'ORFs)
*/


process run_phiSpy {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.run_name}/phispy", mode: 'copy'
    
    container = 'drantoine/phispy:1.0.2'
    // Indispensable car le conteneur ne peut utiliser que les commandes phiSpy
    // Sans ça, il ne reconnaît pas les commandes Bash (mkdir, mv) 
    // et renvoie une erreur "-c eval".
    containerOptions '--entrypoint ""'

    // Gestion de l'erreur 30 (contigs trop courts)
    // Gestion de l'erreur 41 (trop peu d'ORFs)
    errorStrategy { [30, 41].contains(task.exitStatus) ? 'ignore' : 'terminate' }

    input:
    tuple val(sample_id), path(gbk)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: phi_results

    script:
    """
    # Le '|| status=\$?' capture le code d'erreur sans arrêter le script immédiatement
    # En cas d'erreur 30 ou 41, le script devrait normalement écarter les échantillons
    # de l'analyse
    # Dans notre cas, nous souhaitons avant l'arrêt du script créer un fichier
    # README_phispy_skipped.txt qui indiquera que Phispy n'a pu faire l'analyse

    PhiSpy.py ${gbk} \
        -o ${sample_id} \
        --output_choice 33 \
        --threads ${params.phispy_threads} \
        -u ${params.phispy_min_contig_size} \
        -w ${params.phispy_window_size} \
        -n ${params.phispy_number} \
        --extra_dna ${params.phispy_extra_dna} \
        --min_repeat_len ${params.phispy_min_repeat_len}|| status=\$?

    # Si status est vide (cas où PhiSpy a réussi), on le met à 0
    status=\${status:-0}

    # Gestion des codes 30 (Trop court) et 41 (Peu d'ORFs) :
    # Ces codes indiquent que PhiSpy a volontairement ignoré l'échantillon
    # Un dossier est tout de même créer même si PhiSpy a échoué

    if [ "\$status" -eq 30 ] || [ "\$status" -eq 41 ]; then
        echo "Capture de l'erreur \$status pour ${sample_id}. Création du dossier de secours."
        mkdir -p ${sample_id}
        
        # On documente la raison de l'absence de résultats dans le README
        echo "Échantillon ignoré par PhiSpy (Code \$status)" > ${sample_id}/README_phispy_skipped.txt
        echo "Raison : Contigs trop courts ou densité d'ORFs insuffisante pour une prédiction fiable." >> ${sample_id}/README_phispy_skipped.txt
        
        # Le exit 0 permet de mettre le status à 0 si c'était une erreur 30 ou 41
        # Cela permet à Nextflow de considérer le job comme 'SUCCESS' et de publier le dossier
        exit 0
    elif [ "\$status" -ne 0 ]; then
        # Si c'est un autre code d'erreur (ex: crash mémoire), on laisse le pipeline s'arrêter
        exit \$status
    fi
    """
}

// -------------------------------
// DefenseFinder
// -------------------------------

/*
DESCRIPTION
DefenseFinder recherche les systèmes de défense antiviraux (CRISPR, RM, 
Abi, etc.) dans un génome. L'outil prédit les gènes à partir du FASTA avec Prodigal,
puis des profils HMM et MacSyFinder et pour identifier et valider les systèmes.
    --> DefenseFinder va utiliser des HMM afin de détecter
    des features de certaines fonction de défense (Modèle de défense)
    --> MacSyFinder sera appelé par la suite pour
    rechercher les gènes obligatoires puis accesoires
    aux systèmes détectés. Il utilisera également
    la synténie des gènes pour évaluer la distance et l'ordre de ces gènes


SPECIFICITÉ DU CONTAINER ET DES MODÈLES
- Container : 'bjhall/defense-finder:2.0.1-models-2.0.2'
- Modèles : Utilise les modèles CRISPR/Defense 2.0.2.
    --> Modèle Crispr 2.0 : https://github.com/macsy-models/CasFinder/releases/tag/3.1.0

- Option '--skip-model-version-check' : Obligatoire pour utiliser les modèles 
  locaux du container situés dans '/root/.macsyfinder/models'.

PARAMÈTRES RÉGLABLES (via nextflow.config)
- params.defense_workers      : Nombre de cœurs (0 = tous les cœurs, Défaut: 0)
- params.defense_coverage     : Couverture minimale des profils HMM (Défaut: 0.4)
- params.defense_db_type      : Mode d'organisation génomique (Défaut: 'ordered_replicon')
- params.defense_preserve_raw : Conserve les sorties MacSyFinder brutes (Défaut: false)

INPUT
- sample_id : Identifiant de l'échantillon
- fasta     : Séquences nucléotidiques des contigs


OUTPUT 
- Un dossier (${sample_id}/) contenant :
  --> clean.prt : Séquences protéiques que Prodigal a prédit à partir du fasta
  --> clean.prt.idx : Fichier d'indexation pour clean.prt
  --> defense_finder_hmmer.tsv   : Hits trouvés par HMMER
  --> defense_finder_genes.tsv   : Tous les gènes de défense détectés
  --> defense_finder_systems.tsv : Systèmes complets validés par synténie (MacSyFinder)
*/

process run_defenseFinder {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.run_name}/defense", mode: 'copy'

    container = 'bjhall/defense-finder:2.0.1-models-2.0.2'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: def_results

    script:
    """
    defense-finder run ${fasta} \
        -o ${sample_id} \
        -w ${params.defense_workers} \
        -c ${params.defense_coverage} \
        --db-type ${params.defense_db_type} \
        --skip-model-version-check \
        --models-dir /root/.macsyfinder/models \
        ${params.defense_preserve_raw ? '--preserve-raw' : ''}
    """
}

// -------------------------------
// AntiDefenseFinder
// -------------------------------

/*
DESCRIPTION (voir DefenseFinder)
AntiDefenseFinder utilise le programme DefenseFinder (HMM + MacSyFinder) 
mais cible spécifiquement les systèmes d'anti-défense (ex: Anti-CRISPR / Acr)
par l'activation de l'option '--antidefensefinder-only' désactivant
la recherche des systèmes de défense classiques pour se concentrer 
uniquement sur les contre-attaques virales.
*/

process run_antiDefenseFinder {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.run_name}/anti_defense", mode: 'copy'

    container = 'bjhall/defense-finder:2.0.1-models-2.0.2'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: anti_results

    script:
    """
    defense-finder run ${fasta} \
        -o ${sample_id} \
        --antidefensefinder-only \
        -w ${params.defense_workers} \
        -c ${params.defense_coverage} \
        --db-type ${params.defense_db_type} \
        --skip-model-version-check \
        --models-dir /root/.macsyfinder/models \
        ${params.defense_preserve_raw ? '--preserve-raw' : ''}
    """
}

// -------------------------------
// CrisprCasFinder
// -------------------------------

/*
DESCRIPTION
CRISPRCasFinder permet la détection et la classification des systèmes CRISPR-Cas. 
Sa force réside dans sa capacité à identifier les deux composantes du système :
1- Les CRISPR arrays : Successions de répétitions (DR) et de spacer (via Vmatch)
2- Les gènes Cas : Protéines associées identifiées par MacSyFinder et HMMER.
L'outil attribue un "Evidence Level" (de 1 = faible à 4 = excellent) 
aux arrays détectés.

--> Identification des répétitions : Recherche de motifs répétés caractéristiques 
    de l'array CRISPR.
--> Analyse du voisinage (-vi) : Recherche de gènes Cas autour des arrays détectés.
--> Intégration MacSyFinder : Utilise des modèles spécifiques pour classer le 
    système (ex: Type I-E, Type II-A).

PARAMÈTRES RÉGLABLES (via nextflow.config)
- params.crispr_vicinity        : Fenêtre de recherche (bp) autour des arrays pour 
                                  trouver les gènes Cas (Défaut: 20000)
- params.crispr_cpuMacSyFinder  : Nombre de cœurs alloués à MacSyFinder
- params.crispr_keepAll         : Conserve les fichiers temporaires (-keep)

INPUT
- sample_id : Identifiant de l'échantillon
- fasta :Séquences nucléotidiques des contigs. 
    --> Chemin absolu via 'readlink -f', car Perl 
    change souvent de répertoire de travail en interne

OUTPUT 
- Un dossier (${sample_id}/) contenant une arborescence :
  --> result.json            : Toutes les données (positions, types, scores) au format JSON
  --> TSV/                   : Dossier contenant les tableaux récapitulatifs des spacers, 
                               des DR et des gènes Cas
  --> GFF/                   : Annotations avec un fichier par contig
*/

process run_crisprCasFinder {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.run_name}/crispr", mode: 'copy'

    container 'unlhcc/crisprcasfinder:4.2.20'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: cris_results

    script:
    
    """
    mkdir -p crispr

    # -cf et -soFile : Chemins internes au container pointant vers 
    # la configuration et le moteur de recherche Vmatch.

    CRISPRCasFinder.pl \
        -in \$(readlink -f $fasta) \
        -outdir $sample_id \
        -cf /opt/CRISPRCasFinder/config.properties \
        -soFile /opt/CRISPRCasFinder/sel392v2.so \
        -vi ${params.crispr_vicinity} \
        -cpuM ${params.crispr_cpuMacSyFinder} \
        -cas
    """
}


// -------------------------------
// DGR
// -------------------------------

/*
DESCRIPTION
DGR-package (DGRpy) identifie les DGRs. Ces systèmes permettent aux phages et bactéries
de faire muter de manière ciblée certaines protéines (ex : récepteur phage) 
pour s'adapter à leur cible. 

    --> Prédiction d'ORF via Prodigal (intégré)
    --> Identification des RT (Transcriptases Inverses) via HMMER
    --> Alignement TR-VR : Recherche de répétitions imparfaites des TR/VR
    avec des mutations des adénines (séparés par centaines de bases)


OUTPUT (dossier ${sample_id}/)
- ..._DGR_csv.csv        : Résumé global (Score, Mismatches, types de mutations)
- ..._RT_sequences.fasta : Séquences des protéines Reverse Transcriptases identifiées
- ..._hmmer_table.txt    : Statistiques de détection des domaines RT
*/

process run_DGR {

    tag "${sample_id}"

    publishDir "${params.outdir}/${params.run_name}/dgr", mode: 'copy'

    container 'dgr:1.0'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}/"), emit: dgr_results

    script:
    // On définit le nom local pour éviter les problèmes de chemins absolus
    """
# Dans le process prepare data, les fichiers ont été nettoyés
# Tous les fichiers fasta s'appellent dont clean.fasta
# Ceci pose problème à DGR lors de la création de dossiers, 
# le if permet de gérer ce conflit de noms similaires entre les échantillons
    if [ "${fasta}" != "${sample_id}.fasta" ]; then
        cp ${fasta} ${sample_id}.fasta
    fi

    python - <<END
from DGR_package import run_module
import os

fasta_input = "${sample_id}.fasta"

# Exécution du module
run_module.run_file_path(fasta_input)
END

# Par défaut DGR rajoute en suffixe du dossier : _file_directory
# Le mv va éliminer ce suffixe pour simplifier le nom
# Sinon pour un échantillon, si aucun résultat n'est trouvé
# Un dossier sera quand même créée, attendu par le tuple output

    if [ -d "${sample_id}_file_directory" ]; then
        mv ${sample_id}_file_directory ${sample_id}
    else
        mkdir ${sample_id}
    fi
    """
}

process generate_Summary {
   
    publishDir "${params.outdir}/${params.run_name}/_Final_Summary", mode: 'copy'

    container 'biocontainers/pandas:1.5.1_cv1'

    /*
    Le stage as va créer des dossiers temporaires permettant 
    de structurer l'ensemble des fichiers de sortie 
    avec un dossier par outil et par échantillon 
    Pour la suite des étapes, cela facilitera pour Python
    la recherche des fichiers de résultats pour chaque outil
    afin de générer un Summary 
    */


    input:
    val official_ids //Liste des sample_id
    path 'mge_in/*',      stageAs: 'mge_results/*' //Regroupe les sorties de la suite MGEfinder (inferseq) dans 'mge_results'
    path 'int_in/*',      stageAs: 'int_results/*' //Regroupe les sorties d'IntegronFinder dans 'int_results'
    path 'phi_in/*',      stageAs: 'phi_results/*'
    path 'def_in/*',      stageAs: 'def_results/*'
    path 'anti_in/*',     stageAs: 'anti_results/*'
    path 'cris_in/*',     stageAs: 'cris_results/*'
    path 'dgr_in/*',      stageAs: 'dgr_results/*'

    output:
    path "BactPhage-Annot_summary.csv"

    /*
    Transforme la liste d'échantillon ['Sp_1', 'Sp_2', 'Sp_3'] 
    en .txt séparant les échantillons par une \n tel que : 
    Sp_1
    Sp_2
    Sp_3

    Ceci permet de gérer les conflits syntaxiques entre Nextflow et Python
    */

    script:
    def ids_str = official_ids.join('\n') 
    """
    echo "${ids_str}" > ids_list.txt

    python3 <<CODE
import os, glob
import pandas as pd

with open("ids_list.txt") as f:
    ids = [l.strip() for l in f if l.strip()]

# Colonnes qui vont structurer le CSV
metrics = ["MGE", "Integrons", "Defense", "AntiDefense", "Prophages", "CRISPR", "DGR_RT_candidate", "DGR_Complete"]

# Création d'un dictionnaire "vide" pour chaque sample_id
# Si un échantillon n'a aucun résultat, il apparaîtra avec des 0
results = {sid: {m: 0 for m in metrics} for sid in ids}


# Recherche les sample_id dans les différents chemins de chaque outil
# ex : dgr_results/Sp_1/fichier.csv ; dgr_results/Sp_32/fichier.csv
# Il va noter que le premier chemin correspond au sample_id "Sp_1"
# Et que le second correspond au sample_id "Sp_32"
def find_sid(path):
    for sid in ids:
        if f"/{sid}/" in path or path.endswith(f"/{sid}"):
            return sid
    return None

# ------ Lecture des résultats ------

#-------------------------------
# MGE
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers .inferseq.tsv dans 'mge_results'
for f in glob.glob("mge_results/**/*.inferseq.tsv", recursive=True):
    # Fonction décrite précédemment pour savoir à quel échantillon appartient ce fichier
    sid = find_sid(f)

    # len compte le nombre de lignes donc d'éléments détectés et stocke le nombre dans un dictionnaire
    if sid:
        try: results[sid]['MGE'] = len(pd.read_csv(f, sep='\\t'))
        except: pass

#-------------------------------
# INTEGRONS
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers .summary
for f in glob.glob("int_results/**/*.summary", recursive=True):

    sid = find_sid(f)
    # On ignore les lignes commençant par '#' correspondant à l'en-tête
    # Les colonnes 'CALIN', 'complete', 'In0' sont définies
    # Le contenu de chaque colonne est additionné pour le Summary
    if sid:
        try:
            df = pd.read_csv(f, sep='\\s+', comment='#')
            cols = [c for c in ['CALIN', 'complete', 'In0'] if c in df.columns]
            results[sid]['Integrons'] = int(df[cols].sum().sum()) if cols else 0
        except: pass

#-------------------------------
# DEFENSE FINDER
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers .systems.tsv
for f in glob.glob("def_results/**/*_systems.tsv", recursive=True):
    sid = find_sid(f)
    if sid:
        try: results[sid]['Defense'] = len(pd.read_csv(f, sep='\\t'))
        except: pass

#-------------------------------
# ANTIDEFENSE FINDER
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers .systems.tsv
for f in glob.glob("anti_results/**/*_systems.tsv", recursive=True):
    sid = find_sid(f)
    if sid:
        try: results[sid]['AntiDefense'] = len(pd.read_csv(f, sep='\\t'))
        except: pass

#-------------------------------
# PHISPY
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers prophage_coordinates.tsv
for f in glob.glob("phi_results/**/prophage_coordinates.tsv", recursive=True):
    sid = find_sid(f)
    if sid:
        try: results[sid]['Prophages'] = len(pd.read_csv(f, sep='\\t'))
        except: pass

#-------------------------------
# CRISPR
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers Crisprs_REPORT.tsv
for f in glob.glob("cris_results/**/Crisprs_REPORT.tsv", recursive=True):
    sid = find_sid(f)
    if sid:
        try:
            # On utilise skip_blank_lines=True pour éviter que les lignes vides
            # ne soient comptées comme des résultats
            df = pd.read_csv(f, sep='\t', skip_blank_lines=True)
            df = df.dropna(how='all')

            # Chaque ligne = 1 locus CRISPR détecté
            results[sid]['CRISPR'] = len(df)
            
        except:
            # En cas de fichier introuvable ou illisible
            results[sid]['CRISPR'] = 0

#-------------------------------
# DGR
#-------------------------------
# glob.glob : Recherche récursivement tous les fichiers _RT_coordinates_table.csv
# Ce fichier contient uniquement les coordonnées des RT prédites
for f in glob.glob("dgr_results/**/*_RT_coordinates_table.csv", recursive=True):
    sid = find_sid(f)
    if sid:
        try:
            df_rt = pd.read_csv(f)
            results[sid]['DGR_RT_candidate'] = len(df_rt)

            # Si une RT a été trouvée alors on va chercher s'il y'a un système complet dans Sp_1_DGR_csv.csv
            # Au lieu de refaire un for, le programme va partir du principe que
            # Si un fichier Sp_1_RT_coordinates_table.csv a été trouvé dans le dossier Sp_1_genome_output
            # alors il y'aura forcément un fichier Sp_1_DGR_csv.csv dans le dossier au-dessus
            # Le premier replace va tronquer "_genome_output/" pour remonter d'un dossier
            # Le second va remplacer dans le chemin "Sp_1_RT_coordinates_table.csv" par "Sp_1_DGR_csv.csv"
            # Ceci permet de ne pas réutiliser un for et de ne pas reparcourir l'arborescence
            f_complete = f.replace("_genome_output/", "/").replace("_RT_coordinates_table.csv", "_DGR_csv.csv")
            
            # Si le chemin de fichier déduit existe bien
            if os.path.exists(f_complete):
                try:
                    df_sys = pd.read_csv(f_complete)
                    # Si uniquement l'en-tête, len(df_sys) renverra 0, s'il contient des lignes il comptera
                    results[sid]['DGR_Complete'] = len(df_sys)

                # Tolérance de l'erreur si fichier vide
                except pd.errors.EmptyDataError:
                    results[sid]['DGR_Complete'] = 0
                
                # Tolérance de l'erreur si fichier corrompu ou autre cela mettra 0
                except Exception as e:
                    print(f"ALERTE : Erreur inattendue sur {sid} ({f_complete}) : {e}")
                    results[sid]['DGR_Complete'] = 0
            else:
                # Si le fichier f_complete n'existe pas du tout
                print(f"INFO : {sid} - Fichier DGR_csv.csv absent.")
                results[sid]['DGR_Complete'] = 0
        except:
            pass

# ------ Sauvegarde ------
# Création DataFrame à partir du dictionnaire, les clés sample_id deviennent des lignes
df = pd.DataFrame.from_dict(results, orient='index')
df.index.name = 'sample_id'
df[metrics].sort_index().to_csv("BactPhage-Annot_summary.csv")
CODE
    """
}

process generate_GFF {
    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_name}/_GFF_Summary", mode: 'copy'

    input:
    // Ajout de 'def_in' pour les TSV et 'prokka_gff' pour le dictionnaire de positions
    tuple val(sample_id), path('dgr_in'), path('int_in'), path('phi_in'), path('cris_in'), path('def_in'), path('antidef_in'), path('prokka_gff')

    output:
    path "${sample_id}_MGE_map.gff"

    script:
    """
    python3 - <<CODE
import pandas as pd
import glob
import os
import re

all_features = []

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

# ------ Indexation du GFF Prokka ------

prokka_map = {}
if os.path.exists("${prokka_gff}"):
    contig_cds = {}
    with open("${prokka_gff}", 'r') as f:
        for line in f:
            if line.startswith('##FASTA'): break
            if line.startswith('#') or not line.strip(): continue
            p = line.split('\\t')
            if len(p) < 9: continue
            if p[2] == 'CDS':
                contig = p[0]
                if contig not in contig_cds:
                    contig_cds[contig] = []
                # On stocke le début, la fin et le nom du contig
                contig_cds[contig].append((int(p[3]), int(p[4]), p[0]))

    # Pour chaque contig, on numérote les gènes 1, 2, 3... 
    # pour créer des clés comme "C1_1", "C1_2", etc.
    for contig, list_cds in contig_cds.items():
        # On trie par position de début sur le contig
        list_cds.sort() 
        for i, coords in enumerate(list_cds, 1):
            key = f"{contig}_{i}"
            prokka_map[key] = coords

#-------------------------------
# INTEGRONS
#-------------------------------
if os.path.isdir("int_in"):
    for f in glob.glob("int_in/**/*.integrons", recursive=True):
        try:
            df = pd.read_csv(f, sep='\\t', comment='#')
            for _, r in df.iterrows():
                strand = '+' if str(r['strand']) == '1' else '-'
                feat = r['annotation'] if r['annotation'] != 'protein' else 'CDS'
                attr = f"ID={r['element']};Parent={r['ID_integron']};Status={r['type']}"
                all_features.append([r['ID_replicon'], 'IntegronFinder', feat, int(r['pos_beg']), int(r['pos_end']), '.', strand, '.', attr])
        except Exception: pass

#-------------------------------
# DGR
#-------------------------------
if os.path.isdir("dgr_in"):
    for f in glob.glob("dgr_in/**/*_RT_coordinates_table.csv", recursive=True):
        try:
            df = pd.read_csv(f)
            for i, r in df.iterrows():
                strand = '+' if str(r['Direction']) == '1' else '-'
                attr = f"ID=DGR_{i};Clade={r['Clade']};Note=Potential_DGR_reverse_transcriptase"
                all_features.append([r['ID'], 'DGR_Scanner', 'DGR_RT', int(r['Start_Pos']), int(r['Stop_Pos']), '.', strand, '.', attr])
        except Exception: pass

    for f in glob.glob("dgr_in/**/*_DGR_csv.csv", recursive=True):
        try:
            # On vérifie si le fichier n'est pas vide (plus que l'entête)
            if os.path.getsize(f) > 120:
                df = pd.read_csv(f)
                for i, r in df.iterrows():
                    # Ici 'Sequence' est l'ID du contig
                    # Comme on n'a pas les positions exactes Start/Stop du TR/VR dans ce CSV,
                    # On peut mettre des points '.' ou essayer de les inférer si besoin.
                    attr = f"ID=DGR_system_{i};Score={r['Score']};Note=Validated_TR_VR_alignment"
                    # On l'ajoute comme une information supplémentaire sur le contig
                    all_features.append([r['Sequence'], 'DGR_Scanner', 'DGR_system_match', '.', '.', '.', '.', '.', attr])
        except Exception: pass

#-------------------------------
# PHISPY
#-------------------------------
phispy_files = glob.glob("phi_in/**/*.gff*", recursive=True)
for f in phispy_files:
    try:
        with open(f, 'r') as ph:
            for line in ph:
                if line.strip() and not line.startswith("#"):
                    p = line.strip().split('\\t')
                    if len(p) == 9:
                        all_features.append([p[0], p[1], p[2], int(p[3]), int(p[4]), p[5], p[6], p[7], p[8]])
    except Exception: pass

#-------------------------------
# CRISPRCasFinder
#-------------------------------
if os.path.isdir("cris_in"):
    crispr_reports = glob.glob("cris_in/**/Crisprs_REPORT.tsv", recursive=True)
    for f in crispr_reports:
        try:
            df = pd.read_csv(f, sep='\\t')
            for _, r in df.iterrows():
                strand = '+' if r['CRISPRDirection'] == 'Forward' else '-' if r['CRISPRDirection'] == 'Reverse' else '.'
                attr = f"ID={r['CRISPR_Id']};Evidence_Level={r['Evidence_Level']};Spacers={r['Spacers_Nb']}"
                all_features.append([r['Sequence'], 'CRISPRCasFinder', 'CRISPR_array', int(r['CRISPR_Start']), int(r['CRISPR_End']), '.', strand, '.', attr])
        except Exception: pass

#-------------------------------
#  DEFENSE FINDER
#-------------------------------
if os.path.isdir("def_in"):
    def_systems = glob.glob("def_in/**/*_systems.tsv", recursive=True)
    for f in def_systems:
        try:
            df = pd.read_csv(f, sep='\\t')
            for _, r in df.iterrows():
                b_id = str(r['sys_beg']).strip()
                e_id = str(r['sys_end']).strip()
                
                if b_id in prokka_map and e_id in prokka_map:
                    start_pb = prokka_map[b_id][0]
                    end_pb = prokka_map[e_id][1]
                    contig = prokka_map[b_id][2]
                    attr = f"ID={r['sys_id']};Type={r['type']};Subtype={r['subtype']}"
                    all_features.append([contig, 'DefenseFinder', 'defense_system', start_pb, end_pb, '.', '.', '.', attr])
        except Exception:
            pass

#-------------------------------
# ANTIDEFENSE FINDER
#-------------------------------
if os.path.isdir("antidef_in"):
    def_systems = glob.glob("antidef_in/**/*_systems.tsv", recursive=True)
    for f in def_systems:
        try:
            df = pd.read_csv(f, sep='\\t')
            for _, r in df.iterrows():
                b_id = str(r['sys_beg']).strip()
                e_id = str(r['sys_end']).strip()
                
                if b_id in prokka_map and e_id in prokka_map:
                    start_pb = prokka_map[b_id][0]
                    end_pb = prokka_map[e_id][1]
                    contig = prokka_map[b_id][2]
                    attr = f"ID={r['sys_id']};Type={r['type']};Subtype={r['subtype']}"
                    all_features.append([contig, 'AntiDefenseFinder', 'antidefense_system', start_pb, end_pb, '.', '.', '.', attr])
        except Exception:
            pass

# ------ TRI ------
all_features.sort(key=lambda x: (natural_sort_key(str(x[0])), x[3]))

# ------ EN-TÊTE ------
# On commence par la version, puis on va chercher les tailles de contigs dans le GFF Prokka
headers = [
    "##gff-version 3",
    "## Column Header: seqid | source | type | start | end | score | strand | phase | attributes"
]
if os.path.exists("${prokka_gff}"):
    with open("${prokka_gff}", 'r') as f_in:
        for line in f_in:
            if line.startswith('##sequence-region'):
                headers.append(line.strip())
            # Dès qu'on arrive aux gènes (CDS/gene...), on arrête de lire le header
            elif not line.startswith('#'):
                break

# ------ ÉCRITURE ------
out_gff = "${sample_id}_MGE_map.gff"
with open(out_gff, 'w') as gff:
    # 1. Écriture de l'en-tête dynamique
    for h in headers:
        gff.write(h + "\\n")
    
    # 2. Écriture des features triées
    for feat in all_features:
        gff.write("\\t".join(map(str, feat)) + "\\n")
CODE
    """
}


// -------------------------------
// WORKFLOW
// -------------------------------

workflow {

println "=== Pipeline BactPhage-Annot – workflow démarrage ==="

    // ------ Nettoyage et Harmonisation ------
    // On envoie le flux brut du CSV vers le process de préparation
    prepared_ch = PREPARE_DATA(samples_ch)

    // C'est ici qu'on sépare le flux propre pour chaque outil
    // Le GFF est déclaré car il fait partie du tuple mais
    // il servira seulement plus tard pour le Summary et le GFF
    // Ce mapfile fait le lien entre les identifiants de contigs originaux issus de Prokka 
    // et les noms personnalisés définis dans le pipeline
    data = prepared_ch.ready.multiMap { id, fasta, gff_clean, bam, bai, ref, gbk, map_file ->
        phispy:      tuple(id, gbk)
        crispr:      tuple(id, fasta)
        dgr:         tuple(id, fasta)
        integron:    tuple(id, fasta)
        defenseF:    tuple(id, fasta) //C'est le même tuple pour DefenseFinder et AntiDefenseFinder
        mefinder:    tuple(id, ref, bam, bai)
        for_summary: tuple(id, map_file)
    }

    // ------ Mobile Genetic Elements (MGEfinder) ------
    // On utilise data.mefinder (les fichiers propres)
    if (params.run_mgefinder) {
        find_out  = mge_find(data.mefinder)
        pair_out  = mge_pair(find_out)
        infer_out = mge_inferseq(pair_out)
    } else {
        infer_out = Channel.empty()
    }

    // ------ Integrons ------
    if (params.run_integron) {
        integron_out = run_integronFinder(data.integron)
    } else {
        integron_out = Channel.empty()
    } 

    // ------ Prophages (PhiSpy) ------
    if (params.run_phispy) {
        phispy_out = run_phiSpy(data.phispy)
    } else {
        phispy_out = Channel.empty()
    }

    // ---------Defense Systems & Anti-Defense ------
    if (params.run_defense) {
        defense_out      = run_defenseFinder(data.defenseF)
        anti_defense_out = run_antiDefenseFinder(data.defenseF)
    } else {
        defense_out      = Channel.empty()
        anti_defense_out = Channel.empty()
    }

    // ------ CRISPR-Cas ------
    if (params.run_crispr) {
        crispr_out = run_crisprCasFinder(data.crispr)
    } else {
        crispr_out = Channel.empty()
    }

    // ------ DGR ------
    if (params.run_dgr) {
        dgr_out = run_DGR(data.dgr)
    } else {
        dgr_out = Channel.empty()
    }

    // ------ Génération du Résumé ------
    // On utilise data.for_summary (les IDs issus du flux propre)
    all_mappings_ch = data.for_summary.map{ it[0] }.collect() // On ne garde que les IDs

    generate_Summary(
        all_mappings_ch,                              
        infer_out.map{ it[1] }.collect().ifEmpty([]),
        integron_out.map{ it[1] }.collect().ifEmpty([]),
        phispy_out.map{ it[1] }.collect().ifEmpty([]),
        defense_out.map{ it[1] }.collect().ifEmpty([]),
        anti_defense_out.map{ it[1] }.collect().ifEmpty([]),
        crispr_out.map{ it[1] }.collect().ifEmpty([]),
        dgr_out.map{ it[1] }.collect().ifEmpty([])
    )

    // ------ Génération du GFF ------
    // On récupère l'ID (it[0]) et le clean.gff (it[2])
    prokka_gff_ch = PREPARE_DATA.out.ready.map { it -> [it[0], it[2]] }

    ch_for_gff = dgr_out
        .join(integron_out, remainder: true)
        .join(phispy_out, remainder: true)
        .join(crispr_out, remainder: true)
        .join(defense_out, remainder: true)
        .join(anti_defense_out, remainder: true)
        .join(prokka_gff_ch)
        .map { id, dgr, integron, phispy, crispr, deffinder, antideffinder, prokka ->
            return tuple(id, dgr ?: [], integron ?: [], phispy ?: [], crispr ?: [], deffinder ?: [], antideffinder ?: [], prokka)
        }
    
    generate_GFF(ch_for_gff)
}


