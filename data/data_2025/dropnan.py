import pandas as pd

# List of input files
input_files = ["LibraryDesign.csv", "LibraryNatural.csv"]

for file in input_files:
    # Read the CSV
    df = pd.read_csv(file)

    # Drop rows where [RE_norm] is NaN
    df_clean = df.dropna(subset=["RE_norm"])

    # Save cleaned version
    output_file = file.replace(".csv", "_cleaned.csv")
    df_clean.to_csv(output_file, index=False)

    print(f"Cleaned '{file}' -> '{output_file}' with {len(df) - len(df_clean)} rows dropped.")

