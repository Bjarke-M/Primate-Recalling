#!/usr/bin/env python3
"""
For each species, loops over intersection BED files in
  {BASE_DIR}
which are named like:
  {species}_batch{N}_fploidy{F}_mploidy{M}_{unique|nonunique|noalign}.bed
 
Computes total length per chromosome per bed and writes a table with columns:
  SPECIES, BATCH, FPLOIDY, MPLOIDY, BED_TYPE, CHR, N_CALLABLE
 
Usage: python callable_lengths.py
"""
 
import subprocess
import glob
import re
import pandas as pd
 
REF_SPECIES = ["Allenopithecus_nigroviridis_ssp", "Allochrocebus_lhoesti_ssp", "Allochrocebus_preussi_ssp", "Allochrocebus_solatus_ssp", "Alouatta_belzebul_ssp", 
               "Alouatta_caraya_ssp", "Alouatta_discolor_ssp", "Alouatta_juara_ssp", "Alouatta_macconnelli_ssp", "Alouatta_palliata_ssp", "Alouatta_seniculus_ssp", 
               "Aotus_azarai_ssp", "Aotus_griseimembra_ssp", "Aotus_nancymaae_ssp", "Aotus_trivirgatus_ssp", "Aotus_vociferans_ssp", "Arctocebus_calabarensis_ssp", 
               "Ateles_belzebuth_ssp", "Ateles_chamek_ssp", "Ateles_geoffroyi_ssp", "Ateles_marginatus_ssp", "Ateles_paniscus_ssp", "Avahi_laniger_ssp",
               "Avahi_occidentalis_ssp", "Avahi_peyrierasi_ssp", "Brachyteles_hypoxanthus_ssp", "Cacajao_ayresi_ssp", "Cacajao_calvus_ssp", "Cacajao_hosomi_ssp", 
               "Cacajao_melanocephalus_ssp", "Callimico_goeldii_ssp", "Callithrix_geoffroyi_ssp", "Callithrix_jacchus_ssp", "Callithrix_kuhlii_ssp", "Carlito_syrichta_ssp",
               "Cebuella_niveiventris_ssp", "Cebuella_pygmaea_ssp", "Cebus_albifrons_ssp", "Cebus_imitator_ssp", "Cebus_olivaceus_ssp", "Cebus_unicolor_ssp", 
               "Cephalopachus_bancanus_ssp", "Cercocebus_atys_ssp", "Cercocebus_chrysogaster_ssp", "Cercocebus_lunulatus_ssp", "Cercocebus_torquatus_ssp", 
               "Cercopithecus_ascanius_ssp", "Cercopithecus_campbelli_ssp", "Cercopithecus_cephus_ssp", "Cercopithecus_denti_ssp", "Cercopithecus_diana_ssp", 
               "Cercopithecus_hamlyni_ssp", "Cercopithecus_mitis_ssp", "Cercopithecus_mona_ssp", "Cercopithecus_neglectus_ssp", "Cercopithecus_nictitans_ssp", 
               "Cercopithecus_petaurista_ssp", "Cercopithecus_pogonias_ssp", "Cercopithecus_roloway_ssp", "Cercopithecus_wolfi_ssp", "Cheirogaleus_crossleyi_ssp", 
               "Cheirogaleus_major_ssp", "Cheirogaleus_medius_ssp", "Cheirogaleus_sibreei_ssp", "Cheracebus_lucifer_ssp", "Cheracebus_lugens_ssp", "Cheracebus_regulus_ssp", 
               "Cheracebus_torquatus_ssp", "Chiropotes_albinasus_ssp", "Chiropotes_israelita_ssp", "Chiropotes_sagulatus_ssp", "Chlorocebus_aethiops_ssp", 
               "Chlorocebus_cynosuros_ssp", "Chlorocebus_dryas_ssp", "Chlorocebus_pygerythrus_ssp", "Chlorocebus_sabaeus_ssp", "Chlorocebus_tantalus_ssp", 
               "Colobus_angolensis_ssp", "Colobus_guereza_ssp", "Colobus_polykomos_ssp", "Daubentonia_madagascariensis_ssp", "Erythrocebus_patas_ssp", "Eulemur_albifrons_ssp", 
               "Eulemur_cinereiceps_ssp", "Eulemur_collaris_ssp", "Eulemur_coronatus_ssp", "Eulemur_flavifrons_ssp", "Eulemur_fulvus_ssp", "Eulemur_macaco_ssp", 
               "Eulemur_mongoz_ssp", "Eulemur_rubriventer_ssp", "Eulemur_rufifrons_ssp", "Eulemur_rufus_ssp", "Eulemur_sanfordi_ssp", "Galago_moholi_ssp", 
               "Galago_senegalensis_ssp", "Galagoides_demidoff_ssp", "Gorilla_beringei_ssp", "Gorilla_gorilla_ssp", "Hapalemur_aureus_ssp", "Hapalemur_griseus_ssp", 
               "Hapalemur_meridionalis_ssp", "Hapalemur_occidentalis_ssp", "Hoolock_hoolock_ssp", "Hoolock_leuconedys_ssp", "Hylobates_abbotti_ssp", "Hylobates_agilis_ssp", 
               "Hylobates_klossii_ssp", "Hylobates_lar_ssp", "Hylobates_moloch_ssp", "Hylobates_muelleri_ssp", "Hylobates_pileatus_ssp", "Indri_indri_ssp", 
               "Lagothrix_lagotricha_ssp", "Lemur_catta_ssp", "Leontocebus_fuscicollis_ssp", "Leontocebus_nigricollis_ssp", "Leontopithecus_chrysomelas_ssp", 
               "Leontopithecus_rosalia_ssp", "Lepilemur_ankaranensis_ssp", "Lepilemur_dorsalis_ssp", "Lepilemur_mittermeieri_ssp", "Lepilemur_mustelinus_ssp", 
               "Lepilemur_ruficaudatus_ssp", "Lepilemur_sahamalazensis_ssp", "Lepilemur_septentrionalis_ssp", "Lophocebus_aterrimus_ssp", "Loris_lydekkerianus_ssp", 
               "Loris_tardigradus_ssp", "Macaca_arctoides_ssp", "Macaca_assamensis_ssp", "Macaca_brunnescens_ssp", "Macaca_cyclopis_ssp", "Macaca_fascicularis_ssp", 
               "Macaca_fuscata_ssp", "Macaca_hecki_ssp", "Macaca_leonina_ssp", "Macaca_leucogenys_ssp", "Macaca_maura_ssp", "Macaca_mulatta_ssp", "Macaca_nemestrina_ssp", 
               "Macaca_nigra_ssp", "Macaca_nigrescens_ssp", "Macaca_radiata_ssp", "Macaca_siberu_ssp", "Macaca_silenus_ssp", "Macaca_sinica_ssp", "Macaca_sylvanus_ssp", 
               "Macaca_thibetana_ssp", "Macaca_tonkeana_ssp", "Mandrillus_leucophaeus_ssp", "Mandrillus_sphinx_ssp", "Mico_argentatus_ssp", "Mico_humeralifer_ssp", 
               "Mico_humilis_ssp", "Microcebus_arnholdi_ssp", "Microcebus_griseorufus_ssp", "Microcebus_jonahi_ssp", "Microcebus_mittermeieri_ssp", "Microcebus_murinus_ssp", 
               "Microcebus_ravelobensis_ssp", "Microcebus_tavaratra_ssp", "Miopithecus_ogouensis_ssp", "Miopithecus_talapoin_ssp", "Mirza_coquereli_ssp", "Mirza_zaza_ssp", 
               "Nasalis_larvatus_ssp", "Nomascus_annamensis_ssp", "Nomascus_concolor_ssp", "Nomascus_gabriellae_ssp", "Nomascus_leucogenys_ssp", "Nomascus_siki_ssp", 
               "Nycticebus_bengalensis_ssp", "Nycticebus_coucang_ssp", "Otolemur_crassicaudatus_ssp", "Otolemur_garnettii_ssp", "Pan_paniscus_ssp", "Pan_troglodytes_ssp", 
               "Papio_anubis_ssp", "Papio_cynocephalus_ssp", "Papio_hamadryas_ssp", "Papio_kindae_ssp", "Papio_papio_ssp", "Papio_ursinus_ssp", "Perodicticus_potto_ssp", 
               "Piliocolobus_badius_ssp", "Piliocolobus_gordonorum_ssp", "Piliocolobus_kirkii_ssp", "Piliocolobus_tephrosceles_ssp", "Pithecia_albicans_ssp", 
               "Pithecia_chrysocephala_ssp", "Pithecia_hirsuta_ssp", "Pithecia_mittermeieri_ssp", "Pithecia_pissinatti_ssp", "Pithecia_pithecia_ssp", "Pithecia_vanzolinii_ssp", 
               "Plecturocebus_bernhardi_ssp", "Plecturocebus_brunneus_ssp", "Plecturocebus_caligatus_ssp", "Plecturocebus_cinerascens_ssp", "Plecturocebus_cupreus_ssp", 
               "Plecturocebus_donacophilus_ssp", "Plecturocebus_dubius_ssp", "Plecturocebus_grovesi_ssp", "Plecturocebus_hoffmannsi_ssp", "Plecturocebus_miltoni_ssp", 
               "Plecturocebus_moloch_ssp", "Pongo_abelii_ssp", "Pongo_pygmaeus_ssp", "Pongo_tapanuliensis_ssp", "Presbytis_comata_ssp", "Presbytis_melalophos_ssp", 
               "Prolemur_simus_ssp", "Propithecus_coquereli_ssp", "Propithecus_deckenii_ssp", "Propithecus_diadema_ssp", "Propithecus_edwardsi_ssp", "Propithecus_perrieri_ssp", 
               "Propithecus_tattersalli_ssp", "Propithecus_verreauxi_ssp", "Pygathrix_cinerea_ssp", "Pygathrix_nemaeus_ssp", "Pygathrix_nigripes_ssp", 
               "Rhinopithecus_avunculus_ssp", "Rhinopithecus_bieti_ssp", "Rhinopithecus_brelichi_ssp", "Rhinopithecus_roxellana_ssp", "Rhinopithecus_strykeri_ssp", 
               "Saguinus_bicolor_ssp", "Saguinus_geoffroyi_ssp", "Saguinus_imperator_ssp", "Saguinus_inustus_ssp", "Saguinus_labiatus_ssp", "Saguinus_midas_ssp", 
               "Saguinus_mystax_ssp", "Saguinus_oedipus_ssp", "Saimiri_boliviensis_ssp", "Saimiri_cassiquiarensis_ssp", "Saimiri_macrodon_ssp", "Saimiri_oerstedii_ssp", 
               "Saimiri_sciureus_ssp", "Saimiri_ustus_ssp", "Sapajus_apella_ssp", "Semnopithecus_entellus_ssp", "Semnopithecus_hypoleucos_ssp", "Semnopithecus_priam_ssp", 
               "Semnopithecus_schistaceus_ssp", "Symphalangus_syndactylus_ssp", "Tarsius_dentatus_ssp", "Tarsius_lariang_ssp", "Tarsius_wallacei_ssp", "Theropithecus_gelada_ssp", 
               "Trachypithecus_auratus_ssp", "Trachypithecus_cristatus_ssp", "Trachypithecus_ebenus_ssp", "Trachypithecus_francoisi_ssp", "Trachypithecus_geei_ssp", 
               "Trachypithecus_germaini_ssp", "Trachypithecus_hatinhensis_ssp", "Trachypithecus_johnii_ssp", "Trachypithecus_laotum_ssp", "Trachypithecus_obscurus_ssp", 
               "Trachypithecus_phayrei_ssp", "Trachypithecus_pileatus_ssp", "Trachypithecus_poliocephalus_ssp", "Trachypithecus_vetulus_ssp", "Varecia_rubra_ssp", 
               "Varecia_variegata_ssp", "Xanthonycticebus_pygmaeus_ssp"]

BASE_DIR = ""

 
# Filename pattern: {species}_batch{N}_fploidy{F}_mploidy{M}_{bed_type}.bed
PATTERN = re.compile(
    r".+_batch(\d+)_fploidy(\d+)_mploidy(\d+)_(unique|nonunique|noalign)\.bed$"
)
 

 
for species in REF_SPECIES:
    print(species)
    rows = []
    intersections_dir = f"{BASE_DIR}/{species}/"
    OUT_FILE  = f"{BASE_DIR}/{species}/callable_lengths.tsv"
    bed_files = sorted(glob.glob(f"{intersections_dir}/{species}_batch*.bed"))
 
    if not bed_files:
        print(f"WARNING: no intersection beds found in {intersections_dir}")
        continue
 
    for bed_path in bed_files:
        fname = bed_path.split("/")[-1]
        m = PATTERN.match(fname)
        if not m:
            print(f"WARNING: skipping unrecognised filename: {fname}")
            continue
 
        batch, fploidy, mploidy, bed_type = m.group(1), m.group(2), m.group(3), m.group(4)
 
        awk_cmd = f"awk '{{sum[$1] += $3-$2}} END {{for (chr in sum) print chr\"\\t\"sum[chr]}}' {bed_path}"
        result = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)
 
        if result.returncode != 0:
            print(f"WARNING: awk failed for {fname}: {result.stderr.strip()}")
            continue
 
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            chrom, length = line.split("\t")
            rows.append({
                "SPECIES":    species,
                "BATCH":      int(batch),
                "FPLOIDY":    int(fploidy),
                "MPLOIDY":    int(mploidy),
                "BED_TYPE":   bed_type,
                "CHR":        chrom,
                "N_CALLABLE": int(length)
            })
 
    df = pd.DataFrame(rows, columns=["SPECIES", "BATCH", "FPLOIDY", "MPLOIDY", "BED_TYPE", "CHR", "N_CALLABLE"])
    df = df.sort_values(["SPECIES", "BATCH", "FPLOIDY", "MPLOIDY", "BED_TYPE", "CHR"]).reset_index(drop=True)
    df.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"Written: {OUT_FILE} ({len(df)} rows)")
    