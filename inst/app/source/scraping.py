import requests
from bs4 import BeautifulSoup
import pandas as pd
import re
import ast
import time
import json # <-- Added for JSON output

def fetch_survival_rates():
    # Start a session
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15"
    })

    search_url = "https://hbcr-survival.ganjoho.jp/search"
    graph_url = "https://hbcr-survival.ganjoho.jp/graph"

    print("Fetching cancer types from the search page...")
    response = session.get(search_url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')

    cancer_select = soup.find('select', id='cancer_types')
    if not cancer_select:
        print("Error: Could not find the cancer type selection list.")
        return

    cancer_options = []
    for option in cancer_select.find_all('option'):
        val = option.get('value')
        if val and "code" in val:
            text = option.text.strip()
            cancer_options.append((val, text))

    print(f"Total {len(cancer_options)} cancer types to process.")
    all_data = []

    for val, name in cancer_options:
        print(f"-> Fetching data for {name}...")

        # Construct POST payload
        # Changed to "男女" (Both genders) and added "全年齢" (All ages) for fallback
        payload = {
            'cancer_types': val,
            'diagnosis_year': '{"year":"2014-2015","elapsed_year":"5"}',
            'gender1': '{"name":"男女","value":0}', # Both genders
            'stages_40': '{"name":"Ⅳ期","value":40}',
            'age_class_AC99': '{"name":"全年齢","value":"AC99"}', # Fallback overall age
            'age_class_AC00': '{"name":"40未満","value":"AC00"}',
            'age_class_AC40': '{"name":"40代","value":"AC40"}',
            'age_class_AC50': '{"name":"50代","value":"AC50"}',
            'age_class_AC60': '{"name":"60代","value":"AC60"}',
            'age_class_AC70': '{"name":"70代","value":"AC70"}',
            'age_class_AC80': '{"name":"80以上","value":"AC80"}',
            'operation1': '{"name":"全体","value":9}'
        }

        # Send POST request to /graph endpoint
        resp = session.post(graph_url, data=payload)
        html_source = resp.text

        # Extract column names and row data
        columns = re.findall(r"ms_data\.addColumn\('number',\s*'(.*?)'\);", html_source)
        rows_match = re.search(r"ms_data\.addRows\(\[(.*?)\]\);", html_source)

        if not columns or len(columns) <= 1 or not rows_match:
            print(f"   [!] {name}: No graph data")
            continue

        rows_str = "[" + rows_match.group(1) + "]"
        try:
            rows = ast.literal_eval(rows_str)
        except Exception as e:
            continue

        if not rows or len(rows) == 0:
            print(f"   [!] {name}: No graph data (N < 30 for all groups)")
            continue

        if columns[0] == 'year':
            df = pd.DataFrame(rows, columns=columns)
        else:
            df = pd.DataFrame(rows, columns=['year'] + columns)

        df_melted = df.melt(id_vars=['year'], var_name='category', value_name='survival_rate')
        df_melted['cancer_type'] = name
        df_melted = df_melted[df_melted['year'].between(1, 5)]

        def split_category(cat_str):
            parts = cat_str.split(':')
            if len(parts) >= 5:
                return pd.Series(parts[:5]) # label, sex, stage, age_group, surgery
            return pd.Series([cat_str, None, None, None, None])

        df_melted[['label', 'sex', 'stage', 'age_group', 'surgery']] = df_melted['category'].apply(split_category)
        
        all_data.append(df_melted)
        time.sleep(1) # Be gentle with the server

    if not all_data:
        print("\nNo valid data could be retrieved.")
        return

    # Combine data
    final_df = pd.concat(all_data, ignore_index=True)
    final_df = final_df[['cancer_type', 'sex', 'age_group', 'stage', 'surgery', 'year', 'survival_rate']]
    final_df.sort_values(by=['cancer_type', 'age_group', 'year'], inplace=True)

    # Output CSV for reference
    final_df.to_csv("survival_rates_stage4_2014-2015_tidy.csv", index=False, encoding="utf-8-sig")
    print("\n[Done] CSV saved as: survival_rates_stage4_2014-2015_tidy.csv")

    # =========================================================
    # Generate JSON format directly for R Shiny integration
    # =========================================================
    output_json = {}
    for cancer_name, group_cancer in final_df.groupby('cancer_type'):
        output_json[cancer_name] = {}
        # Group by age_group (e.g., '全年齢', '40代', '50代'...)
        for age, group_age in group_cancer.groupby('age_group'):
            rates = group_age.sort_values('year')['survival_rate'].tolist()
            if len(rates) == 5:
                output_json[cancer_name][age] = rates
                
    with open("Data_age_survival_5_year.json", "w", encoding="utf-8") as f:
        json.dump(output_json, f, ensure_ascii=False, indent=2)
    print("[Done] JSON data saved as: Data_age_survival_5_year.json")

if __name__ == "__main__":
    fetch_survival_rates()

