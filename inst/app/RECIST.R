classification_rules <- list(
  "Blood and lymphatic system disorders" = c(
    "anemia", "neutropenia", "thrombocytopenia", "leukocytosis", "eosinophilia",
    "neutrophil count decreased", "platelet count decreased", "white blood cell decreased",
    "lymphocyte count decreased", "bone marrow hypocellular", "hemolysis",
    "disseminated intravascular coagulation", "hemolytic uremic syndrome",
    "myelodysplastic syndrome", "febrile neutropenia", "cd4 lymphocytes decreased"
  ),
  
  "Cardiac disorders" = c(
    "heart failure", "myocarditis", "myocardial infarction", "atrial fibrillation",
    "atrial flutter", "left ventricular systolic dysfunction", "pericarditis",
    "pericardial effusion", "pericardial tamponade", "cardiac arrest",
    "sinus bradycardia", "supraventricular tachycardia", "ventricular arrhythmia",
    "paroxysmal atrial tachycardia", "restrictive cardiomyopathy", "conduction disorder",
    "endocarditis infective", "ejection fraction decreased", "chest pain - cardiac",
    "palpitations"
  ),
  
  "Gastrointestinal disorders" = c(
    "diarrhea", "nausea", "vomiting", "constipation", "abdominal pain", "colitis",
    "enterocolitis", "gastritis", "esophagitis", "pancreatitis", "mucositis oral",
    "pharyngeal mucositis", "colonic perforation", "jejunal perforation", "rectal perforation",
    "small intestinal perforation", "duodenal perforation", "gastric perforation",
    "ileal perforation", "lower gastrointestinal hemorrhage", "upper gastrointestinal hemorrhage",
    "cecal hemorrhage", "duodenal hemorrhage", "gastric hemorrhage", "colonic hemorrhage",
    "ileal hemorrhage", "jejunal hemorrhage", "pancreatic hemorrhage", "pharyngeal hemorrhage",
    "small intestinal obstruction", "ileal obstruction", "colonic obstruction",
    "duodenal obstruction", "rectal obstruction", "duodenal ulcer", "gastric ulcer",
    "colonic ulcer", "rectal ulcer", "anal ulcer", "stomach pain", "pharyngeal fistula",
    "rectal fistula", "anal fistula", "esophageal fistula", "enterovesical fistula",
    "pancreatic fistula", "small intestinal stenosis", "duodenal stenosis", "rectal stenosis",
    "gastric stenosis", "enterocolitis infectious", "appendicitis", "cholecystitis",
    "ileus", "ascites", "chylous ascites", "dysphagia", "dyspepsia", "dysgeusia",
    "hiccups", "esophageal varices hemorrhage", "small intestinal mucositis", "laryngeal mucositis",
    "tracheal mucositis", "bile duct stenosis", "intra-abdominal hemorrhage", "abdominal distension",
    "rectal pain", "anal pain", "oral pain", "anorexia", "appetite loss", "fecal incontinence",
    "colonic fistula"
  ),
  
  "Renal and urinary disorders" = c(
    "proteinuria", "acute kidney injury", "chronic kidney disease", "urinary tract infection",
    "creatinine increased", "hematuria", "nephrotic syndrome", "kidney infection",
    "renal hemorrhage", "cystitis noninfective", "urinary tract obstruction",
    "glucosuria", "urethral infection"
  ),
  
  "Respiratory, thoracic and mediastinal disorders" = c(
    "pneumonitis", "dyspnea", "pleural effusion", "bronchopulmonary hemorrhage",
    "lung infection", "pulmonary fibrosis", "pneumothorax", "adult respiratory distress syndrome",
    "hypoxia", "wheezing", "respiratory failure", "pulmonary edema", "bronchial stricture",
    "pleural infection", "atelectasis", "bronchial infection", "tracheal hemorrhage",
    "productive cough", "interstitial pneumonia", "tracheitis", "bronchopleural fistula",
    "pulmonary hypertension", "pulmonary embolism", "interstitial lung disease", 
    "cough", "hoarseness", "aspiration", "nasal congestion", "sinusitis", 
    "upper respiratory infection", "mediastinal infection", "vital capacity abnormal",
    "pharyngeal necrosis"
  ),
  
  "Infections and infestations" = c(
    "sepsis", "gallbladder infection", "lung infection", "bone infection", "mucosal infection",
    "soft tissue infection", "bacteremia", "enterocolitis infectious", "biliary tract infection",
    "hepatic infection", "kidney infection", "gum infection", "catheter related infection",
    "abdominal infection", "mediastinal infection", "small intestine infection",
    "lymph gland infection", "salivary gland infection", "pelvic infection", "pancreas infection",
    "infective myositis", "pleural infection", "encephalitis infection", "bronchial infection",
    "tooth infection", "peritoneal infection", "anorectal infection", "thrush",
    "biliary tract infection", "breast infection", "urethral infection",
    "esophageal infection", "splenic infection", "cranial nerve infection",
    "corneal infection", "nail infection", "skin infection", "pharyngitis",
    "meningitis", "shingles"
  ),
  
  "Skin and subcutaneous tissue disorders" = c(
    "eczema", "rash acneiform", "erythema multiforme", "stevens-johnson syndrome",
    "papulopustular rash", "urticaria", "palmar-plantar erythrodysesthesia syndrome",
    "rash maculo-papular", "paronychia", "erythroderma", "rash", "skin ulceration",
    "bullous dermatitis", "dry skin", "photosensitivity", "nail changes", "alopecia",
    "pruritus", "rash pustular", "cheilitis", "skin hyperpigmentation", "nail discoloration",
    "nail loss", "hair color changes", "pain of skin", "skin induration", "rash acnelform",
    "toxic epidermal necrolysis", "dermatitis radiation", "purpura"
  ),
  
  "Nervous system disorders" = c(
    "peripheral sensory neuropathy", "peripheral motor neuropathy", "encephalopathy",
    "neuralgia", "myasthenia gravis", "gait disturbance", "seizure", "dysesthesia",
    "acoustic nerve disorder nos", "guillain-barre syndrome", "ischemia cerebrovascular",
    "stroke", "dysarthria", "glossopharyngeal nerve disorder", "akathisia",
    "extrapyramidal disorder", "recurrent laryngeal nerve palsy", "headache",
    "cognitive disturbance", "optic nerve disorder", "brachial plexopathy",
    "depressed level of consciousness", "vertigo", "facial nerve disorder",
    "reversible posterior leukoencephalopathy syndrome", "vagus nerve disorder",
    "arachnoiditis", "oculomotor nerve disorder", "dizziness", "confusion",
    "leukoencephalopathy", "ataxia", "generalized muscle weakness", "syncope",
    "intracranial hemorrhage"
  ),
  
  "Endocrine disorders" = c(
    "hypophysitis", "adrenal insufficiency", "hypothyroidism", "hypopituitarism",
    "hyperthyroidism", "hyperparathyroidism", "hypoparathyroidism"
  ),
  
  "Metabolism and nutrition disorders" = c(
    "hyponatremia", "hypokalemia", "hyperglycemia", "hyperkalemia", "hypomagnesemia",
    "dehydration", "glucose intolerance", "weight loss", "hypocalcemia", "hypotension",
    "hypertriglyceridemia", "acidosis", "hyperphosphatemia", "hypophosphatemia",
    "hypoglycemia", "hypoalbuminemia", "hypercalcemia", "weight gain", "cholesterol high"
  ),
  
  "General disorders and administration site conditions" = c(
    "malaise", "edema limbs", "allergic reaction", "infusion related reaction",
    "fatigue", "back pain", "vasovagal reaction", "cytokine release syndrome",
    "epistaxis", "pain", "localized edema", "generalized edema", "tumor hemorrhage",
    "multi-organ failure", "tumor lysis syndrome", "disease progression", "lymph node pain",
    "fall", "injection site reaction", "wound complication", "edema face",
    "chest wall pain", "edema trunk", "flu like symptoms", "bone pain",
    "hot flashes", "non-cardiac chest pain", "periorbital edema", "death nos",
    "infusion site extravasation", "fever"
  ),
  
  "Musculoskeletal and connective tissue disorders" = c(
    "myalgia", "arthritis", "myositis", "osteonecrosis of jaw", "arthralgia",
    "osteonecrosis", "muscle cramp", "infective myositis", "fracture", "rhabdomyolysis",
    "ankle fracture", "muscle weakness lower limb", "toothache"
  ),
  
  "Psychiatric disorders" = c(
    "depression", "anxiety", "insomnia", "hallucinations"
  ),
  
  "Vascular disorders" = c(
    "hypertension", "thromboembolic event", "thrombotic thrombocytopenic purpura",
    "vasculitis", "portal vein thrombosis", "arterial thromboembolism",
    "peripheral ischemia", "portal hypertension", "superior vena cava syndrome",
    "lymphedema"
  ),
  
  "Eye disorders" = c(
    "retinal detachment", "keratitis", "blurred vision", "retinopathy", "uveitis",
    "eyelid function disorder", "vision decreased", "corneal ulcer", "watering eyes",
    "conjunctivitis"
  ),
  
  "Ear and labyrinth disorders" = c(
    "hearing impaired"
  ),
  
  "Hepatobiliary disorders" = c(
    "aspartate aminotransferase increased", "alanine aminotransferase increased",
    "serum amylase increased", "lipase increased", "blood bilirubin increased",
    "alkaline phosphatase increased", "alt", "ggt increased", "alt elevation",
    "alt increased", "elevated alt", "hepatic failure", "hepatic necrosis",
    "perforation bile duct", "gallbladder perforation", "hepatic pain",
    "hepatic hemorrhage", "hepatobiliary disorder", "transeaminase elevation",
    "sinusoidal obstruction syndrome"
  ),
  
  "Investigations" = c(
    "cpk increased", "blood lactate dehydrogenase increased", "electrocardiogram qt corrected interval prolonged",
    "blood antidiuretic hormone abnormal"
  ),
  
  "Immune system disorders" = c(
    "anaphylaxis", "autoimmune disorder"
  ),
  
  "Reproductive system and breast disorders" = c(
    "vaginal fistula", "uterine hemorrhage", "uterine fistula", "vaginal inflammation"
  ),
  
  "Neoplasms benign, malignant and unspecified" = c(
    "leukemia secondary to oncology chemotherapy", "treatment related secondary malignancy"
  ),
  
  "Injury, poisoning and procedural complications" = c(
    "alcohol intolerance"
  )
)