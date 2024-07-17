from coreg import CoReg
import pandas as pd
from sklearn.preprocessing import LabelEncoder


if __name__ == "__main__":
    print("run_coreg running")
    k1 = 3
    k2 = 5
    p1 = 2
    p2 = 7
    max_iters = 250
    pool_size = 1000

    cr = CoReg(k1, k2, p1, p2, max_iters, pool_size, verbose=True)

    # reading in and preprocessing data
    df = pd.read_csv('../merged_data_sub_10.csv')
    columns_to_drop = ["SDMVPSU", "SDMVSTRA", "WTMEC4YR", "SDDSRVYR", "RIDSTATR", "WTSCY4YR"]
    df = df.drop(columns=columns_to_drop)
    # Encode categorical variables
    label_encoders = {}
    for column in ["RIAGENDR", "RIDRETH2", "RIDEXPRG"]:
        le = LabelEncoder()
        df[column] = le.fit_transform(df[column].astype(str))
        label_encoders[column] = le
        
    # separating and printing shapes
    labeled_data = df[df['simulated_mGFR'].notna()]
    unlabeled_data = df[df['simulated_mGFR'].isna()]
    unlabeled_data = unlabeled_data.drop(columns=['simulated_mGFR'])
    



    cr.set_datasets(labeled_data, unlabeled_data)

    cr.train()