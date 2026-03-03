file_clinical = "/Users/ikegami/Desktop/length_bias/20250710_池上先生ダミーデータ/clinical_data_whole.qs"
file_variant = "/Users/ikegami/Desktop/length_bias/20250710_池上先生ダミーデータ/variant_data_whole.qs"
file_drug = "/Users/ikegami/Desktop/length_bias/20250710_池上先生ダミーデータ/drug_data_whole.qs"
data_clinical = QS_READ(8,file_clinical)
data_variant = QS_READ(8,file_variant)
data_drug = QS_READ(8,file_drug)

UIDs = unique(data_clinical$C.CAT調査結果.基本項目.ハッシュID)
sample_no = 10000
downsampled_IDs = sample(UIDs, sample_no)

data_clinical_downsampled = data_clinical %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% downsampled_IDs)
data_variant_downsampled = data_variant %>% dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% downsampled_IDs)
data_drug_downsampled = data_drug %>% dplyr::filter(ID %in% downsampled_IDs)

file_clinical_downsampled = "/Users/ikegami/Desktop/length_bias/20250710_池上先生ダミーデータ/clinical_data_whole_downsampled.qs"
file_variant_downsampled = "/Users/ikegami/Desktop/length_bias/20250710_池上先生ダミーデータ/variant_data_whole_downsampled.qs"
file_drug_downsampled = "/Users/ikegami/Desktop/length_bias/20250710_池上先生ダミーデータ/drug_data_whole_downsampled.qs"

QS_SAVE(8,data_clinical_downsampled,file_clinical_downsampled)
QS_SAVE(8,data_variant_downsampled,file_variant_downsampled)
QS_SAVE(8,data_drug_downsampled,file_drug_downsampled)
