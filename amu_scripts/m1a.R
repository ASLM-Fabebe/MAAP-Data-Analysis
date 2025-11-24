#loading antibiotic data resources

antibiotic_classes_amc <- c(
  "Aminoglycosides",
  "Amphenicols",
  "Antibiotics - Intestinal antiinfectives",
  "Carbapenems",
  "Combinations of drugs for treatment of tuberculosis",
  "Combinations of sulfonamides and trimethoprim, incl. derivatives",
  "First-generation cephalosporins",
  "Fluoroquinolones",
  "Fourth-generation cephalosporins",
  "Glycopeptide antibacterials",
  "Hydrazides",
  "Imidazole derivatives",
  "Lincosamides",
  "Macrolides",
  "Nitrofuran-derivatives",
  "Nitroimidazole derivatives",
  "Other agents against amoebiasis and other protozoal diseases",
  "Other drugs for treatment of tuberculosis",
  "Other cephalosporins",
  "Oxazolidinones",
  "Penicillins",
  "Polymyxins",
  "Rifamycins",
  "Second-generation cephalosporins",
  "Sulfonamides",
  "Tetracyclines",
  "Third-generation cephalosporins",
  "Triazole and tetrazole derivatives",
  "Tetracyclines",
  "Penicillins with extended spectrum",
  "Beta-lactamase sensitive penicillins",
  "Beta-lactamase resistant penicillins",
  "Beta-lactamase inhibitors",
  "Penicillins, incl. beta-lactamase inhibitors",
  "Monobactams",
  "Other cephalosporins and penems",
  "Trimethoprim and derivatives",
  "Short-acting sulfonamides",
  "Intermediate-acting sulfonamides",
  "Long-acting sulfonamides",
  "Streptogramins",
  "Streptomycins",
  "Other aminoglycosides",
  "Other quinolones",
  "Combinations of antibacterials",
  "Steroid antibacterials",
  "Nitrofuran derivatives",
  "Other antibacterials"
)

# Exclusion lists
inhibitors <- c(
  "Clavulanic acid", "Sulbactam", "Tazobactam", "Avibactam",'clavulanic',
  "Vaborbactam", "Relebactam", "Zidebactam", "Nacubactam",
  "PAβN", "D13-9001", "EDTA", "cilastatin", "betamipron", 'inhibitor'
)

other_meds_components <- c('benzathine',  #this is a stabilizer
                           'procaine') #this is an anesthetic)

other_non_antibiotics <- c(
  "diloxanide", "praziquantel", "antimalarial", "ors",'antibiotic','antifungal','cephalosporins', 'sodium',
  "i.v", "tabs", "cap", "syrup", "acid", "other",'oral', 'pack','plus', 'tablet', 'tablets', 'capsule','capsules'
)

combination_molecules <- c("chlortetracycline,demeclocycline,tetracycline",
                           "pivampicillin,pivmecillinam",
                           "benzathinebenzylpenicillin,benzylpenicillin,procainebenzylpenicillin",
                           "ampicillin,cloxacillin",
                           "ampicillin,oxacillin",
                           "ampicillin,flucloxacillin",
                           "amoxicillin,cloxacillin",
                           "sulfacarbamide,sulfadiazine,sulfadimidine",
                           "sulfamethoxazole,trimethoprim",
                           "sulfadiazine,trimethoprim",
                           "sulfametrole,trimethoprim",
                           "sulfadiazin,tetroxoprim",
                           "sulfamerazin,trimethoprim",
                           "metronidazole,spiramycin",
                           "levofloxacin,ornidazole",
                           "azithromycin,fluconazole,secnidazole",
                           "ofloxacin,ornidazole",
                           "ciprofloxacin,metronidazole",
                           "ciprofloxacin,tinidazole",
                           "ciprofloxacin,ornidazole",
                           "norfloxacin,tinidazole",
                           "azithromycin,cefixime",
                           "furazolidone,metronidazole",
                           "diloxanide,metronidazole",
                           "isoniazid,rifampicin",
                           "isoniazid,pyrazinamide,rifampicin",
                           "ethambutol,isoniazid,pyrazinamide,rifampicin",
                           "ethambutol,isoniazid,rifampicin")  ##to handle combination products

# Combine exclusions
exclude_extra <- tolower(c(inhibitors, other_non_antibiotics))

#who ref
ddd_ref <- read_excel('amc_resources/ab_molecules_amc.xlsx') %>%
  mutate(name_route=tolower(name_route)) %>%
  distinct(name_route, .keep_all = T)

#atc molecules
atcs_cleaned <- read_excel('amc_resources/ab_molecules.xlsx') %>%
  mutate(original_entry=tolower(trimws(original_entry))) %>%
  distinct(original_entry, .keep_all = T)

# Clean ab dictionary
antibiotics_mol_dict <- trimws(tolower(atcs_cleaned$original_entry))

#strength units reference
units_ref <- read_excel('amc_resources/strength_units_reference.xlsx')

