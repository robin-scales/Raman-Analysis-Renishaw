import pandas as pd


def load(filepath: str):
    # Initial load of csv
    df = pd.read_csv(filepath, sep='\t+', header=0, engine='python') # , comment='#'
    # Get the header list
    res = df.columns.tolist()
    # Enumerate through the headers
    for i, col in enumerate(res):
        # If the header name contains a hash
        if "#" in col:
            # Remove the hash
            res[i] = col.replace("#", "")
    # Rname the columns
    df.columns = res
    print(df)
    print(res)

    match res:
        case ['Wave', 'Intensity']:
            print("Loaded in single point")
            scanType = "SP"
        case ['X', 'Y', 'Wave', 'Intensity']:
            print("Loaded in UNprocessed peak location map")
            scanType = "uMAP"
        case ['XList', 'YList', 'IList']:
            print("Loaded in processed peak location map")
            scanType = "pMAP"
        case _:
          print("UNKNOWN CASE")
          scanType = "MISSING"
    
    return df, scanType
