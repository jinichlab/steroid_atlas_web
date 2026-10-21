"""Main-organism breakdown of the atlas — from the organism strings we already have.

No UniProt lineage fetch needed. Every entry in proteins.csv already has an
'organism' field like "Sus scrofa (Pig)" or "Aspergillus fumigatus".
We extract the genus (first word), classify the top-N genera into broad
taxa via a hand-curated dictionary (~200 common genera covering most of
the atlas), and produce a 3-panel figure:

  A. Superkingdom / kingdom (Eukaryota-broad / Fungi / Bacteria / Other)
  B. Broad lineage (Mammals, Fish, Birds, Reptiles/Amphibians, Invertebrates,
     Fungi, Bacteria, Plants, Protists, Other/unclassified)
  C. Top-15 species

Reads:  ../data/proteins.csv
Writes: ../analysis/organism_summary.txt         human-readable numbers
        ../analysis/organism_figure.png / .pdf   3-panel figure
"""
from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
OUT_TXT = output("organism_summary.txt")
OUT_PNG = output("organism_figure.png")
OUT_PDF = output("organism_figure.pdf")


# Hand-curated genus → broad-lineage mapping.  Covers common lab / model /
# livestock organisms plus the top microbial and invertebrate genera.  Anything
# not in the dictionary lands in "Other / unclassified".
GENUS_TO_TAXON = {
    # === Mammals ===
    **{g: "Mammals" for g in (
        "Leptonychotes Lobodon Hydrurga Ommatophoca "
        "Phodopus Neotoma Onychomys Peromyscus Sigmodon "
        "Oryctolagus Lepus Ochotona Sylvilagus "
        "Sciurus Urocitellus Ictidomys Marmota Spermophilus Cynomys Xerus Tamias "
        "Callospermophilus Ammospermophilus "
        "Castor Aplodontia Pedetes "
        "Jaculus Dipodomys Chaetodipus Perognathus "
        "Microtus Arvicola Myodes Neodon Ondatra Ellobius Lasiopodomys "
        "Catagonus Tayassu Pecari Phacochoerus Hylochoerus "
        "Moschus Odocoileus Rangifer Cervus Alces Dama Muntiacus Capreolus Elaphurus "
        "Phyllostomus Molossus Hipposideros Miniopterus Corynorhinus Lasiurus Nyctalus "
        "Vespertilio Eptesicus Tadarida "
        "Theropithecus Nycticebus Cheirogaleus Varecia Indri "
        "Homo Pan Gorilla Pongo Nomascus Hylobates Symphalangus "
        "Macaca Papio Chlorocebus Mandrillus Cercocebus Piliocolobus Colobus "
        "Rhinopithecus Trachypithecus Nasalis Presbytis Semnopithecus "
        "Callithrix Cebus Aotus Saimiri Ateles Alouatta Sapajus Cebuella "
        "Microcebus Otolemur Nycticebus Propithecus Prolemur Eulemur Lemur Daubentonia "
        "Carlito Tarsius Tupaia "
        "Mus Rattus Cricetulus Mesocricetus Cavia Cricetus Peromyscus Neotoma "
        "Meriones Marmota Ictidomys Spermophilus Xerus Nannospalax Heterocephalus "
        "Fukomys Chinchilla Octodon Dolichotis Hydrochoerus Erethizon Myocastor "
        "Sus Bos Bison Ovis Capra Rupicapra Ovibos Bubalus Syncerus "
        "Camelus Vicugna Lama Sus Odocoileus Rangifer Cervus Alces Dama "
        "Capreolus Muntiacus Elaphurus Antilocapra Giraffa Okapia "
        "Equus Ceratotherium Diceros Tapirus "
        "Canis Vulpes Lycaon Otocyon Nyctereutes Speothos "
        "Ursus Ailuropoda Ailurus Panda Melursus Tremarctos "
        "Felis Panthera Puma Prionailurus Neofelis Lynx Acinonyx Leopardus "
        "Mustela Enhydra Lutra Meles Taxidea Gulo Martes Neovison Ictonyx Mellivora "
        "Suricata Herpestes Cryptoprocta Fossa Nandinia Genetta Paradoxurus "
        "Odobenus Zalophus Otaria Arctocephalus Callorhinus Neophoca Eumetopias "
        "Halichoerus Phoca Pusa Erignathus Cystophora Mirounga Monachus Neomonachus "
        "Balaena Balaenoptera Megaptera Eschrichtius Physeter Kogia Ziphius "
        "Berardius Hyperoodon Mesoplodon Delphinapterus Monodon Orcinus "
        "Tursiops Delphinus Stenella Lagenorhynchus Neophocaena Phocoena "
        "Lipotes Inia Sotalia Pontoporia "
        "Manis Loxodonta Elephas Mammuthus Trichechus Dugong "
        "Erinaceus Sorex Talpa Condylura Scalopus "
        "Rhinolophus Myotis Pipistrellus Miniopterus Pteropus Rousettus Artibeus "
        "Desmodus Eptesicus Nyctalus Vespertilio Corynorhinus Lasiurus "
        "Choloepus Bradypus Dasypus Chaetophractus Tolypeutes "
        "Phascolarctos Trichosurus Vombatus Macropus Notamacropus Wallabia "
        "Petaurus Sarcophilus Dasyurus Antechinus Perameles Isoodon "
        "Monodelphis Didelphis Marmosa "
        "Ornithorhynchus Tachyglossus Zaglossus "
        "Chrysochloris Elephantulus Orycteropus Procavia Loxodonta "
    ).split()},

    # === Fish (Actinopterygii, Chondrichthyes, Cyclostomata) ===
    **{g: "Fish" for g in (
        "Labrus Choerodon Halichoeres Thalassoma Xyrichtys Coris Symphodus Notolabrus "
        "Gadus Melanogrammus Merlangius Pollachius Trisopterus Merluccius Micromesistius Boreogadus "
        "Dissostichus Pagothenia Pogonophryne Gymnodraco Chionodraco Champsocephalus Trematomus "
        "Notothenia Aethotaxis Cryodraco "
        "Albula Elops Megalops "
        "Mola Ranzania Masturus "
        "Acanthochromis "
        "Scomber Rastrelliger Auxis Katsuwonus Thunnus Sarda Scomberomorus "
        "Zoarces Lycodes Melanostigma Lycenchelys "
        "Cirrhinus Labeo Osteochilus Cirrhinus "
        "Parambassis Ambassis Chanda "
        "Austrofundulus Millerichthys Rachovia Terranatos Simpsonichthys "
        "Paramormyrops Mormyrus Gnathonemus Petrocephalus Marcusenius Campylomormyrus "
        "Mugilogobius Rhinogobius Chaenogobius Odontobutis Rhyacichthys "
        "Phoxinus Rhodeus Acheilognathus Tanakia Rhinichthys Notemigonus "
        "Knipowitschia Pomatoschistus Aphia "
        "Scleropages Osteoglossum Arapaima Pantodon "
        "Myripristis Holocentrus Sargocentron Neoniphon "
        "Mastacembelus Macrognathus "
        "Cyclopterus "
        "Esox Novumbra Umbra Dallia "
        "Hucho Brachymystax Parahucho "
        "Anabas Anabantidae Betta Trichogaster Trichopodus Osphronemus Colisa "
        "Sphaeramia Ostorhinchus Cheilodipterus Apogon "
        "Seriola Trachurus Caranx Selene Coryphaena Rachycentron "
        "Salarias Blennius Chasmodes "
        "Electrophorus Gymnotus Sternopygus Apteronotus "
        "Gouania Lepadogaster "
        "Echeneis "
        "Neogobius Gobius Gobiodon Neogobius Ponticola Ctenogobius "
        "Channa Ophicephalus "
        "Triplophysa Barbatula Cobitis Misgurnus "
        "Oncorhynchus Salmo Salvelinus Coregonus Thymallus "
        "Cyprinus Sinocyclocheilus Danio Carassius Ctenopharyngodon "
        "Puntigrus Barbus Labeo Rhinichthys Gobio Rutilus "
        "Oryzias Poecilia Xiphophorus Gambusia Fundulus Nothobranchius "
        "Kryptolebias Aphyosemion Cyprinodon Menidia Atherinops "
        "Oreochromis Sarotherodon Astatotilapia Neolamprologus Haplochromis "
        "Maylandia Pundamilia Metriaclima Amphilophus Cichlasoma "
        "Pterophyllum Symphysodon Geophagus Uaru "
        "Astyanax Colossoma Piaractus Serrasalmus Pygocentrus Metynnis Silurana "
        "Larimichthys Miichthys Nibea Larimichthys Sciaenops Argyrosomus "
        "Perca Sander Etheostoma Percina "
        "Morone Micropterus Lepomis "
        "Scophthalmus Paralichthys Solea Cynoglossus Verasper "
        "Takifugu Tetraodon Fugu Diodon Sphoeroides "
        "Latimeria Callorhinchus Rhinochimaera "
        "Petromyzon Lampetra Ichthyomyzon "
        "Anguilla Conger Ophichthus Muraena "
        "Amphiprion Dascyllus Chromis Pomacentrus Stegastes "
        "Notothenia Champsocephalus Trematomus Chaenocephalus Chionodraco "
        "Boleophthalmus Periophthalmus Monopterus "
        "Hippocampus Syngnathus Nerophis "
        "Cottoperca Cottus Gasterosteus Pungitius Culaea "
        "Ictalurus Silurus Clarias Pangasius Hemibagrus Bagre Ameiurus Mystus "
        "Osmerus Mallotus Salangichthys "
        "Denticeps Alosa Clupea Engraulis Sardina Sardinella "
        "Erpetoichthys Polypterus Chanos Lepisosteus Amia Acipenser Huso "
        "Aldrovandia Notacanthus "
        "Lates Oplegnathus Sparus Dicentrarchus "
        "Boleophthalmus Periophthalmus Monopterus "
        "Anabarilius Rhincodon Carcharodon Sphyrna Prionace Squalus "
    ).split()},

    # === Birds ===
    **{g: "Birds" for g in (
        "Malurus Ptilonorhynchus Sericulus Menura "
        "Cairina Netta Marmaronetta Tadorna Chenonetta Dendrocygna "
        "Sylvia Curruca Phylloscopus Acrocephalus Locustella Cettia Hippolais "
        "Nyctibius Steatornis Podargus Caprimulgus Chordeiles Aegotheles "
        "Hirundo Riparia Delichon Progne Petrochelidon Stelgidopteryx "
        "Centropus Crotophaga Guira Coccyzus "
        "Loxia Carduelis Spinus Carpodacus Pyrrhula Chloris Coccothraustes "
        "Chloropsis Aegithina Chlorocichla "
        "Lepidothrix Manacus Chiroxiphia Pipra Corapipo "
        "Prunella Anthus Motacilla Dendronanthus "
        "Crypturellus Nothocercus Tinamus Rhynchotus Eudromia Nothoprocta "
        "Vidua Estrilda Amandava Lonchura Uraeginthus "
        "Zosterops Taeniopygia Poephila Nothoprocta Strigops Geospiza Camarhynchus "
        "Cactornis Certhidea Platyspiza "
        "Sturnella Icterus Molothrus Agelaius Xanthocephalus Dolichonyx Quiscalus "
        "Mimus Toxostoma Dumetella "
        "Nipponia Threskiornis Plegadis "
        "Gallus Meleagris Coturnix Numida Phasianus Bonasa Tetrao Alectoris "
        "Anas Anser Cygnus Branta Aix Aythya Somateria Bucephala Mergus "
        "Balearica Grus Gruella Antigone "
        "Aptenodytes Pygoscelis Eudyptes Eudyptula Spheniscus "
        "Struthio Casuarius Dromaius Rhea Apteryx Tinamus Nothura "
        "Passer Fringilla Serinus Emberiza Zonotrichia Junco Melospiza Spizella "
        "Corvus Pica Nucifraga Perisoreus Cyanocitta "
        "Sturnus Turdus Erithacus Ficedula Muscicapa Phoenicurus Saxicola "
        "Bombycilla Regulus Certhia Sitta Parus Poecile Cyanistes "
        "Falco Aquila Buteo Circus Accipiter Milvus Haliaeetus "
        "Larus Sterna Fulmarus Diomedea Puffinus Oceanodroma Pelecanus Fregata "
        "Egretta Ardea Nycticorax Bubulcus Ciconia Mycteria "
        "Rallus Fulica Porphyrio Gallinula "
        "Charadrius Vanellus Pluvialis Numenius Tringa Actitis Calidris Limosa Scolopax "
        "Columba Zenaida Streptopelia "
        "Psittacus Amazona Ara Cacatua Nestor Nymphicus Melopsittacus "
        "Tyto Athene Strix Bubo Ninox Otus Asio "
        "Alcedo Ceryle Merops Coracias Upupa Bucorvus Buceros "
        "Ramphastos Aulacorhynchus Selenidera "
        "Picus Dryocopus Colaptes Dendrocopos Sphyrapicus "
        "Muscicapa Ficedula Sialia Turdus Catharus Hylocichla "
        "Troglodytes Certhia Sitta Poecile Baeolophus Aegithalos "
    ).split()},

    # === Reptiles ===
    **{g: "Reptiles" for g in (
        "Laticauda Enhydrina Aipysurus Emydocephalus Pelamis "
        "Gopherus Chelonoidis Testudo Geochelone Aldabrachelys Manouria "
        "Pseudonaja Notechis Oxyuranus Acanthophis Pseudechis Dendroaspis Hydrophis "
        "Salvator Tupinambis Callopistes Dracaena "
        "Alligator Crocodylus Gavialis Caiman Osteolaemus Mecistops Tomistoma "
        "Chelonia Caretta Eretmochelys Dermochelys Trachemys Chrysemys Pseudemys "
        "Terrapene Emys Malaclemys Chelydra "
        "Pelodiscus Apalone Trionyx Pelusios Podocnemis Pelomedusa "
        "Podarcis Lacerta Zootoca Timon Anguis Elgaria Ophisaurus "
        "Varanus Heloderma "
        "Sceloporus Anolis Basiliscus Iguana Amblyrhynchus Conolophus "
        "Hemitheconyx Eublepharis Gekko Hemidactylus Tarentola Phelsuma "
        "Naja Bungarus Ophiophagus Micrurus Elapsoidea "
        "Python Pythonidae Boa Corallus Eunectes Morelia "
        "Crotalus Sistrurus Agkistrodon Bothrops Bitis Vipera Cerastes Echis "
        "Thamnophis Nerodia Coluber Elaphe Pantherophis Lampropeltis "
        "Sphenodon "
    ).split()},

    # === Amphibians ===
    **{g: "Amphibians" for g in (
        "Geotrypetes Microcaecilia Rhinatrema Ichthyophis Caecilia Typhlonectes "
        "Engystomops Physalaemus Pseudopaludicola Leptodactylus Odontophrynus Ceratophrys "
        "Eleutherodactylus Craugastor Pristimantis Adelotus Limnodynastes "
        "Pelobates Spea Scaphiopus "
        "Leptobrachium Spea Scaphiopus Alytes Bombina Discoglossus "
        "Xenopus Silurana Bufo Rhinella Rana Lithobates Pelophylax "
        "Hyla Litoria Osteopilus Dendrobates Oophaga Phyllobates Ranitomeya "
        "Ambystoma Notophthalmus Cynops Pleurodeles Triturus Salamandra "
        "Amphiuma Cryptobranchus Andrias Necturus "
    ).split()},

    # === Invertebrates ===
    **{g: "Invertebrates" for g in (
        "Rotaria Adineta Brachionus Philodina "
        "Drosophila Anopheles Aedes Culex Musca Ceratitis Bactrocera Rhagoletis "
        "Nasonia Apis Bombus Megachile Osmia Vespa Polistes Vespula Formica "
        "Camponotus Solenopsis Atta Acromyrmex Pogonomyrmex Wasmannia Cardiocondyla "
        "Bombyx Manduca Spodoptera Helicoverpa Plutella Chilo Ostrinia Lymantria "
        "Tribolium Sitophilus Dendroctonus Ips Anoplophora Callosobruchus Rhynchophorus "
        "Coccinella Harmonia "
        "Sipha Acyrthosiphon Myzus Aphis Pemphigus Rhopalosiphum "
        "Rhodnius Triatoma Cimex "
        "Blattella Blatta Periplaneta Zootermopsis Reticulitermes Coptotermes "
        "Locusta Schistocerca Melanoplus Gryllus Acheta Teleogryllus "
        "Aplysia Lymnaea Biomphalaria Physella Bulinus Helix Cornu "
        "Octopus Sepia Loligo Nautilus Euprymna "
        "Crassostrea Ostrea Mytilus Mya Ruditapes Argopecten Pinctada Haliotis "
        "Homarus Nephrops Palinurus Panulirus Penaeus Litopenaeus Marsupenaeus "
        "Macrobrachium Callinectes Cancer Chionoecetes Portunus Uca "
        "Daphnia Artemia Triops "
        "Ixodes Amblyomma Rhipicephalus Dermacentor Hyalomma Boophilus "
        "Tetranychus Panonychus Sarcoptes Dermatophagoides "
        "Caenorhabditis Pristionchus Ascaris Brugia Wuchereria Onchocerca Loa Trichuris "
        "Ancylostoma Necator Trichinella Strongyloides Haemonchus Meloidogyne Globodera "
        "Schistosoma Fasciola Clonorchis Opisthorchis Echinococcus Taenia Hymenolepis "
        "Nematostella Hydra Aiptasia Acropora Anthopleura Aurelia "
        "Amphimedon Ephydatia Suberites "
        "Strongylocentrotus Lytechinus Paracentrotus Sphaerechinus Arbacia "
        "Branchiostoma Ciona Halocynthia Molgula Botryllus Oikopleura "
    ).split()},

    # === Fungi ===
    **{g: "Fungi" for g in (
        "Aspergillus Penicillium Talaromyplus Talaromyces Monascus "
        "Fusarium Trichoderma Verticillium Metarhizium Beauveria Cordyceps Ophiocordyceps "
        "Botrytis Sclerotinia Monilinia Stemphylium Alternaria Cochliobolus Bipolaris "
        "Cladosporium Curvularia Setosphaeria Zymoseptoria Mycosphaerella Pseudocercospora "
        "Neurospora Sordaria Podospora Chaetomium Thielavia Myceliophthora "
        "Magnaporthe Pyricularia Colletotrichum "
        "Saccharomyces Kluyveromyces Debaryomyces Pichia Ogataea Yarrowia Kazachstania "
        "Zygosaccharomyces Torulaspora Naumovozyma Vanderwaltozyma Lachancea "
        "Candida Meyerozyma Clavispora Yamadazyma "
        "Schizosaccharomyces Neolentinus "
        "Cryptococcus Filobasidium Vishniacozyma Rhodotorula Sporidiobolus "
        "Malassezia Trichosporon "
        "Ustilago Sporisorium Melanotaenium Anthracocystis Moesziomyces "
        "Puccinia Melampsora Uromyces Cronartium Endocronartium "
        "Trichophyton Microsporum Epidermophyton Arthroderma "
        "Histoplasma Paracoccidioides Coccidioides Blastomyces Emmonsia Talaromyces "
        "Rhizopus Mucor Absidia Cunninghamella Lichtheimia Rhizomucor Mortierella "
        "Coprinopsis Coprinus Agaricus Pleurotus Lentinula Flammulina Volvariella "
        "Serpula Coniophora Postia Antrodia Trametes Phanerochaete "
        "Ustilaginoidea Malassezia "
        "Wallemia "
    ).split()},

    # === Plants ===
    **{g: "Plants" for g in (
        "Arabidopsis Thellungiella Boechera Brassica Raphanus Camelina Sinapis Eutrema "
        "Oryza Zea Sorghum Setaria Panicum Miscanthus Hordeum Triticum Aegilops Secale "
        "Avena Brachypodium Lolium Festuca "
        "Glycine Medicago Trifolium Lotus Cicer Pisum Phaseolus Vigna Cajanus Arachis "
        "Solanum Nicotiana Capsicum Petunia Lycopersicon "
        "Manihot Ricinus Jatropha Hevea "
        "Populus Salix Betula Corylus Quercus Castanea Fagus Alnus "
        "Vitis Prunus Malus Pyrus Fragaria Rubus Rosa Sorbus "
        "Citrus Poncirus Fortunella "
        "Coffea Camellia Theobroma Corchorus Gossypium "
        "Cucumis Cucurbita Citrullus Momordica "
        "Amaranthus Beta Spinacia Chenopodium "
        "Selaginella Physcomitrium Marchantia Sphagnum Anthoceros "
        "Chlamydomonas Volvox Chlorella Micromonas "
    ).split()},

    # === Bacteria ===
    **{g: "Bacteria" for g in (
        "Escherichia Salmonella Shigella Klebsiella Enterobacter Citrobacter Cronobacter "
        "Serratia Proteus Providencia Morganella Yersinia Erwinia Pantoea "
        "Bacillus Geobacillus Anoxybacillus Paenibacillus Lysinibacillus Brevibacillus "
        "Lactobacillus Lactiplantibacillus Ligilactobacillus Companilactobacillus "
        "Limosilactobacillus Latilactobacillus Lentilactobacillus Levilactobacillus "
        "Enterococcus Streptococcus Lactococcus Leuconostoc Pediococcus Oenococcus "
        "Weissella Tetragenococcus "
        "Staphylococcus Macrococcus "
        "Listeria Brochothrix "
        "Corynebacterium Rhodococcus Nocardia Mycobacterium Mycolicibacterium Gordonia "
        "Bifidobacterium Actinomyces Actinoplanes Streptomyces Kitasatospora Micromonospora "
        "Amycolatopsis Nocardiopsis Saccharomonospora Saccharothrix "
        "Pseudomonas Acinetobacter Moraxella Psychrobacter Haemophilus Actinobacillus "
        "Pasteurella Mannheimia Aggregatibacter Bordetella Achromobacter "
        "Neisseria Kingella "
        "Chlamydia Chlamydophila Chlamydophila Waddliya "
        "Rickettsia Orientia Ehrlichia Anaplasma Wolbachia "
        "Coxiella Legionella Francisella "
        "Vibrio Aliivibrio Photobacterium Salinivibrio Enterovibrio Grimontia "
        "Helicobacter Campylobacter Wolinella Sulfurospirillum Arcobacter "
        "Burkholderia Ralstonia Cupriavidus Paraburkholderia Pandoraea "
        "Xanthomonas Stenotrophomonas Xylella Pseudoxanthomonas "
        "Rhodobacter Paracoccus Rhodovulum Sphingomonas Sphingopyxis Novosphingobium "
        "Sinorhizobium Bradyrhizobium Rhizobium Agrobacterium Mesorhizobium "
        "Nitrosomonas Nitrobacter Nitrospira Thiobacillus "
        "Thermus Deinococcus Meiothermus "
        "Cyanothece Nostoc Anabaena Prochlorococcus Synechocystis Synechococcus Microcystis "
        "Bacteroides Prevotella Porphyromonas Parabacteroides Alloprevotella Bacteroidetes "
        "Fusobacterium Leptotrichia "
        "Clostridium Clostridioides Ruminococcus Faecalibacterium Roseburia Eubacterium "
        "Blautia Anaerostipes Butyrivibrio Coprococcus Dorea Christensenella Peptoclostridium "
        "Peptoniphilus Peptostreptococcus Anaerococcus Finegoldia Fusicatenibacter "
        "Akkermansia Verrucomicrobium Prosthecobacter "
    ).split()},

    # === Archaea ===
    **{g: "Archaea" for g in (
        "Methanococcus Methanocaldococcus Methanosarcina Methanobrevibacter "
        "Methanothermobacter Methanospirillum Methanoregula Methanocella Methanoculleus "
        "Thermococcus Pyrococcus Palaeococcus Thermotoga "
        "Sulfolobus Saccharolobus Aeropyrum Pyrobaculum Thermoproteus Desulfurococcus "
        "Halobacterium Haloferax Haloarcula Halococcus Halogeometricum Natronomonas "
        "Nitrososphaera Cenarchaeum Nitrosopumilus "
    ).split()},

    # === Protists ===
    **{g: "Protists" for g in (
        "Plasmodium Babesia Theileria Toxoplasma Neospora Sarcocystis Eimeria Cryptosporidium "
        "Trypanosoma Leishmania Crithidia Bodo "
        "Entamoeba Naegleria Acanthamoeba Balamuthia "
        "Giardia Trichomonas Tritrichomonas "
        "Paramecium Tetrahymena Ichthyophthirius Colpoda Stylonychia Oxytricha "
        "Chlamydomonas Volvox Ostreococcus Micromonas Bathycoccus "
        "Phaeodactylum Thalassiosira Nannochloropsis Ectocarpus Aureococcus Fragilariopsis "
        "Phytophthora Pythium Peronospora Hyaloperonospora Plasmopara Bremia "
        "Emiliania Isochrysis Chrysochromulina Prymnesium "
        "Cyanidioschyzon Galdieria Porphyra Chondrus Gracilaria "
        "Dictyostelium Polysphondylium Physarum "
    ).split()},
}


def broad_group(genus: str) -> str:
    if not genus:
        return "Other / unclassified"
    return GENUS_TO_TAXON.get(genus, "Other / unclassified")


def main() -> int:
    print(f"Reading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    p["organism"] = p["organism"].fillna("").astype(str)
    n_total = len(p)
    print(f"  {n_total:,} atlas entries")

    # Extract species (strip the "(common name)" suffix)
    p["species"] = p["organism"].str.split("(").str[0].str.strip()
    p["genus"] = p["organism"].str.split().str[0]
    p["broad_group"] = p["genus"].apply(broad_group)
    n_species = p["species"].nunique()
    n_genera = p["genus"].nunique()
    print(f"  {n_species:,} distinct species · {n_genera:,} distinct genera")

    n_unclassified = int((p["broad_group"] == "Other / unclassified").sum())
    print(f"  {n_unclassified:,} entries unclassified by our genus dictionary "
          f"({n_unclassified/n_total*100:.1f}%)")

    # Summary text
    lines = []
    lines.append("=" * 72)
    lines.append(f"ORGANISM BREAKDOWN — all {n_total:,} atlas proteins")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"Distinct species (organism strings, sans common name): {n_species:,}")
    lines.append(f"Distinct genera (first word of organism):              {n_genera:,}")
    lines.append("")
    lines.append("Broad lineage groups (from manual genus dictionary):")
    for grp, n in p["broad_group"].value_counts().items():
        lines.append(f"  {n:>6,}  ({n/n_total*100:5.1f}%)  {grp}")
    lines.append("")
    lines.append("Top 20 species:")
    for sp, n in p["species"].value_counts().head(20).items():
        lines.append(f"  {n:>5}  {sp}")
    lines.append("")
    lines.append("Top 20 genera:")
    for g, n in p["genus"].value_counts().head(20).items():
        lines.append(f"  {n:>5}  {g}")
    OUT_TXT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_TXT.name}")

    # Figure
    BAR = "#1F4B99"
    OTHER = "#94A3B8"

    fig = plt.figure(figsize=(13.0, 4.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.4])
    axB = fig.add_subplot(gs[0, 0])
    axC = fig.add_subplot(gs[0, 1])

    # PANEL A (was superkingdom — now broad groups since we're not fetching lineage)
    order = ["Mammals", "Fish", "Birds", "Reptiles", "Amphibians", "Invertebrates",
             "Fungi", "Bacteria", "Archaea", "Plants", "Protists", "Other / unclassified"]
    grp = p["broad_group"].value_counts().reindex(order).dropna()
    grp = grp[grp > 0].sort_values(ascending=True)
    colors = [OTHER if g == "Other / unclassified" else BAR for g in grp.index]
    ypos = np.arange(len(grp))
    axB.barh(ypos, grp.values, color=colors, height=0.7,
             edgecolor="white", linewidth=0.4)
    for y, (name, count) in zip(ypos, grp.items()):
        pct = 100 * count / n_total
        axB.text(count * 1.03, y, f" {count:,} ({pct:.1f}%)",
                 va="center", ha="left", fontsize=9, color="#111827")
    axB.set_yticks(ypos)
    axB.set_yticklabels(grp.index, fontsize=9)
    axB.tick_params(axis="y", length=0, pad=4)
    axB.tick_params(axis="x", labelsize=9)
    axB.set_xlim(0, grp.max() * 1.4)
    axB.set_xlabel("Proteins", fontsize=10)
    axB.set_title("A. Broad lineage groups (from genus dictionary)",
                  fontsize=11, loc="left", pad=8)
    axB.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    for s in ("top", "right", "left"):
        axB.spines[s].set_visible(False)
    axB.grid(True, axis="x", alpha=0.25, linewidth=0.5)

    # PANEL C — top-15 species
    top_sp = p["species"].value_counts().head(15)
    ypos = np.arange(len(top_sp))[::-1]
    axC.barh(ypos, top_sp.values, color=BAR, height=0.7,
             edgecolor="white", linewidth=0.4)
    for y, (name, count) in zip(ypos, top_sp.items()):
        pct = 100 * count / n_total
        axC.text(count * 1.02, y, f" {count:,}  ({pct:.2f}%)",
                 va="center", ha="left", fontsize=9, color="#111827")
    axC.set_yticks(ypos)
    axC.set_yticklabels(top_sp.index, fontsize=9)
    axC.tick_params(axis="y", length=0, pad=4)
    axC.tick_params(axis="x", labelsize=9)
    axC.set_xlim(0, top_sp.max() * 1.30)
    axC.set_xlabel("Proteins", fontsize=10)
    axC.set_title("B. Top 15 species", fontsize=11, loc="left", pad=8)
    axC.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    for s in ("top", "right", "left"):
        axC.spines[s].set_visible(False)
    axC.grid(True, axis="x", alpha=0.25, linewidth=0.5)

    fig.suptitle(f"Organism composition of the atlas  "
                 f"(n = {n_total:,} proteins · {n_species:,} species · {n_genera:,} genera)",
                 fontsize=11.5, x=0.01, ha="left", y=1.04)

    for f in (OUT_PNG, OUT_PDF):
        fig.savefig(f, dpi=300, bbox_inches="tight")
        print(f"Wrote {f.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
